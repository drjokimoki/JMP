/*
  RS vs Mean, Time-Series DGP (AR(1) factor + idio) — CRN-correct, optimized C
  + CSV output with bias–variance (rs_results.csv)
  + Iterate all subset sizes z = 2..N
  + RS is UNWEIGHTED here (uniform random subsets — no |t|-weights)
  + STREAMING prints: baseline first, then each z as soon as it finishes

  Build:
    gcc -O3 -march=native -ffast-math -funroll-loops -DUSE_OMP -fopenmp rs_full_z_to_csv.c -o rs_full_z_to_csv -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef USE_OMP
#include <omp.h>
#endif

/* =========================
   Parameters (match your R)
   ========================= */
static const unsigned MAIN_SEED = 778;
static const int n_rep = 500;
static const int T_len_list[] = {500, 100};
static const int T_len_list_len = 2;
static const int N = 50;
static const double sigma_e = 1.0;
static const int K_RS = 1000;
static const double rho_list[] = {0, 0.8};
static const int rho_list_len = 2;
static const double R2_list[] = {0.3, 0.5, 0.1};
static const int R2_list_len = 3;
static const int OOS_LEN_FIXED = 20;      // fixed 20 OOS points
static const double ridge = 1e-10;

/* =========================
   RNG (xoshiro128** + Box–Muller polar)
   ========================= */
typedef struct { unsigned s[4]; } xrng;
typedef struct { int has_spare; double spare; } norm_state;

static unsigned splitmix32(unsigned *x){
    unsigned z = (*x += 0x9E3779B9u);
    z = (z ^ (z >> 16)) * 0x85EBCA6Bu;
    z = (z ^ (z >> 13)) * 0xC2B2AE35u;
    return z ^ (z >> 16);
}
static inline unsigned rotl(const unsigned x, int k){ return (x << k) | (x >> (32 - k)); }
static void xrng_seed(xrng *r, unsigned seed){
    unsigned x = seed ? seed : 0x12345678u;
    r->s[0]=splitmix32(&x); r->s[1]=splitmix32(&x); r->s[2]=splitmix32(&x); r->s[3]=splitmix32(&x);
}
static unsigned xrng_next(xrng *r){
    const unsigned result = rotl(r->s[1] * 5u, 7) * 9u;
    const unsigned t = r->s[1] << 9;
    r->s[2]^=r->s[0]; r->s[3]^=r->s[1]; r->s[1]^=r->s[2]; r->s[0]^=r->s[3];
    r->s[2]^=t; r->s[3]=rotl(r->s[3],11);
    return result;
}
static double xrng_uniform(xrng *r){ return ((xrng_next(r) >> 8)+1.0)/16777217.0; }
static double xrng_normal(xrng *r, norm_state *st){
    if (st->has_spare){ st->has_spare=0; return st->spare; }
    double u,v,s;
    do {
        u = 2.0 * xrng_uniform(r) - 1.0;
        v = 2.0 * xrng_uniform(r) - 1.0;
        s = u*u + v*v;
        if (s >= 1.0 || s == 0.0) continue;
        double mul = sqrt(-2.0 * log(s) / s);
        st->spare = v * mul; st->has_spare = 1;
        return u * mul;
    } while (1);
}

/* =========================
   Small linear algebra
   ========================= */
static int chol_inplace(double *A, int n){
    for (int i=0;i<n;i++){
        for (int j=0;j<=i;j++){
            double sum = A[i*n+j];
            for (int k=0;k<j;k++) sum -= A[i*n+k]*A[j*n+k];
            if (i==j){
                sum += ridge;
                if (sum <= 0.0) return -1;
                A[i*n+j] = sqrt(sum);
            } else {
                A[i*n+j] = sum / A[j*n+j];
            }
        }
        for (int j=i+1;j<n;j++) A[i*n+j] = 0.0;
    }
    return 0;
}
static void fwdsolve(const double *L, double *b, int n){
    for (int i=0;i<n;i++){
        double s=b[i];
        for (int k=0;k<i;k++) s -= L[i*n+k]*b[k];
        b[i] = s / L[i*n+i];
    }
}
static void backsolveT(const double *L, double *y, int n){
    for (int i=n-1;i>=0;i--){
        double s=y[i];
        for (int k=i+1;k<n;k++) s -= L[k*n+i]*y[k];
        y[i] = s / L[i*n+i];
    }
}
static int chol_solve(double *A, double *b, int n){
    if (chol_inplace(A,n)!=0) return -1;
    fwdsolve(A,b,n); backsolveT(A,b,n);
    return 0;
}

/* =========================
   Symmetric-lower helpers
   ========================= */
static inline double sget_lower(const double *S, int n, int i, int j){
    int I = (i>=j)? i : j;
    int J = (i>=j)? j : i;
    return S[I*n + J];
}
static double quadform_lower(const double *S_lower, int n, const double *w){
    double q=0.0;
    for (int i=0;i<n;i++){
        q += w[i]*w[i]*sget_lower(S_lower,n,i,i);
        for (int j=0;j<i;j++){
            double sij = sget_lower(S_lower,n,i,j);
            q += 2.0 * w[i]*w[j]*sij;
        }
    }
    return q;
}

/* =========================
   DGP helpers
   ========================= */
static void toeplitz_cov(double rho, int p, double *Sigma){
    for (int i=0;i<p;i++){
        Sigma[i*p+i]=1.0;
        double v = rho;
        for (int j=i+1;j<p;j++){ Sigma[i*p+j]=Sigma[j*p+i]=v; v*=rho; }
    }
}
static int chol_copy(const double *A, int n, double *L){ memcpy(L,A,sizeof(double)*n*n); return chol_inplace(L,n); }

/* Draw T rows of Z ~ N(0, Σ) using Σ's Cholesky (lower) */
static void mvn_rows(int T, int p, const double *L, xrng *r, norm_state *st, double *out){
    for (int t=0;t<T;t++){
        double z[128];
        for (int j=0;j<p;j++) z[j] = xrng_normal(r, st);
        for (int i=0;i<p;i++){
            double s=0.0;
            for (int k=0;k<=i;k++) s += L[i*p+k]*z[k];
            out[(size_t)t*p + i] = s;
        }
    }
}

/* ---------- Stationary AR(1) scalar with innovation variance s2 ---------- */
static void gen_AR1_scalar_innovvar(int T, double phi, double s2, xrng *r, norm_state *st, double *x){
    double sd = sqrt(s2);
    /* Start from the invariant distribution so the population-R2 calibration
       remains valid throughout the simulated sample. */
    x[0] = sqrt(s2 / (1.0 - phi*phi)) * xrng_normal(r, st);
    for (int t=1;t<T;t++) x[t] = phi*x[t-1] + sd*xrng_normal(r, st);
}

/* ---------- Stationary vector AR(1) with innovation covariance Sigma ---------- */
static void gen_E_AR1_innovcov(int T, int p, double phi,
                               const double *Sigma_chol, xrng *r, norm_state *st, double *E){
    double *Z = (double*)malloc(sizeof(double)*(size_t)T*p);
    mvn_rows(T, p, Sigma_chol, r, st, Z);               // Z_t ~ N(0, Σ)
    const double stationary_scale = 1.0 / sqrt(1.0 - phi*phi);
    for (int j=0;j<p;j++) E[j] = stationary_scale * Z[j]; /* invariant draw */
    for (int t=1;t<T;t++){
        double *e_t = E + (size_t)t*p;
        double *e_m = E + (size_t)(t-1)*p;
        double *z_t = Z + (size_t)t*p;
        for (int j=0;j<p;j++) e_t[j] = phi*e_m[j] + z_t[j];
    }
    free(Z);
}

/* ---------- i.i.d. white noise with variance s2 ---------- */
static void gen_iid_noise(int T, double s2, xrng *r, norm_state *st, double *x){
    double sd = sqrt(s2);
    x[0] = 0.0;
    for (int t=1; t<T; t++) x[t] = sd * xrng_normal(r, st);
}

/* =========================
   RS helpers (Gram + submatrix)
   ========================= */
static void build_grams(const double *X_tr, const double *y_tr, int T_tr, int N,
                        double *S_lower, double *b, double *s1, double *Ty){
    for (int i=0;i<N;i++){ b[i]=0.0; s1[i]=0.0; }
    for (int i=0;i<N;i++) for (int j=0;j<=i;j++) S_lower[i*N+j]=0.0;
    double ty=0.0;
    for (int t=0;t<T_tr;t++){
        const double *xt = X_tr + (size_t)t*N;
        double yt = y_tr[t];
        ty += yt;
        for (int i=0;i<N;i++){
            double xi = xt[i];
            b[i]  += xi * yt;
            s1[i] += xi;
            for (int j=0;j<=i;j++){
                S_lower[i*N+j] += xi * xt[j];
            }
        }
    }
    *Ty = ty;
}
static void gather_lower(const double *S_lower, int N, const int *sub, int z, double *S_sub){
    for (int i=0;i<z;i++){
        int ii = sub[i];
        for (int j=0;j<=i;j++){
            int jj = sub[j];
            int I = ii>=jj ? ii : jj;
            int J = ii>=jj ? jj : ii;
            S_sub[i*z + j] = S_lower[I*N + J];
        }
        for (int j=i+1;j<z;j++) S_sub[i*z + j] = 0.0;
    }
}
static void sample_k_no_replace(xrng *r, int n, int k, int *pool, int *sub){
    for (int i=0;i<n;i++) pool[i]=i;
    for (int i=0;i<k;i++){
        /* The modulo introduces negligible bias here and, unlike the previous
           floating-point mapping, can never produce the out-of-range index n. */
        int j = i + (int)(xrng_next(r) % (unsigned)(n - i));
        int tmp = pool[i]; pool[i] = pool[j]; pool[j] = tmp;
        sub[i] = pool[i];
    }
}

/* =========================
   Σ_u builder and sum-to-one optimal weights
   ========================= */
static void build_sigma_u_lower(const double *X, const double *y_tr,
                                const double *gamma, int T_tr, int N,
                                double *Sigma_u_lower)
{
    for (int i=0;i<N;i++) for (int j=0;j<=i;j++) Sigma_u_lower[i*N+j]=0.0;
    for (int t=0;t<T_tr;t++){
        const double *xt = X + (size_t)t*N;
        double yt = y_tr[t];
        double urow[128];
        for (int i=0;i<N;i++) urow[i] = yt - gamma[i]*xt[i];       /* base error: y - f_i */
        for (int i=0;i<N;i++){
            for (int j=0;j<=i;j++){
                Sigma_u_lower[i*N+j] += urow[i]*urow[j];
            }
        }
    }
    double invT = 1.0 / (double)T_tr;
    for (int i=0;i<N;i++) for (int j=0;j<=i;j++) Sigma_u_lower[i*N+j] *= invT;
}

static int sum_to_one_weights_from_cov_lower(const double *Sigma_lower, int n,
                                             double *w_out, double *Q_out)
{
    double *A = (double*)malloc(sizeof(double)*n*n);
    for (int i=0;i<n;i++){
        for (int j=0;j<=i;j++) A[i*n+j] = Sigma_lower[i*n+j];
        for (int j=i+1;j<n;j++) A[i*n+j] = 0.0;
    }
    double *b = (double*)malloc(sizeof(double)*n);
    for (int i=0;i<n;i++) b[i]=1.0;

    int ok = chol_solve(A, b, n); /* b = Sigma^{-1} 1 */
    if (ok!=0){ free(A); free(b); return -1; }

    double denom=0.0; for (int i=0;i<n;i++) denom += b[i];
    if (denom<=0.0){ free(A); free(b); return -2; }

    for (int i=0;i<n;i++) w_out[i] = b[i] / denom;
    if (Q_out) *Q_out = 1.0 / denom;      /* GMV Q */

    free(A); free(b);
    return 0;
}

/* =========================
   Main
   ========================= */
int main(void){
#ifdef USE_OMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif
    printf("Running with %d thread(s)\n", num_threads);
    setvbuf(stdout, NULL, _IONBF, 0);   /* immediate prints */

    /* open CSV and write header */
    FILE *fp = fopen("rs_results.csv","w");
    if (!fp){ perror("fopen rs_results.csv"); return 1; }
    fprintf(fp,
        "T_len,rho,R2_target,phi_f,phi_eta,phi_u,z,"
        "rmsfe_base,rmsfe_rs,rel_rmsfe,"
        "rmse_pooled_base,rmse_pooled_rs,"
        "bias_base,var_base,bias_rs,var_rs,"
        "Q_ols,Q_rs_bag,Q_rs_subavg,"
        "L_ols,L_rs_bag,L_rs_subavg\n");
    fflush(fp);

    /* phi triplets (no duplicates) */
    const int P=2;
    const double phi_f_arr[2]   = {0.3, 0.7};
    const double phi_eta_arr[2] = {0.7, 0.3};
    const double phi_u_arr[2]   = {0.0, 0.0};  /* reported only */

    /* beta_raw: first 10 ones */
    double beta_raw[N]; for (int i=0;i<N;i++) beta_raw[i] = (i<10)?1.0:0.0;

    typedef struct { int T_len; double rho, R2_target, phi_f, phi_eta, phi_u; } Scenario;

    /* Optional baseline T-scaling mode. The paper-design run remains the
       default; set RSM_T_SCALING=1 for a fixed-DGP experiment varying only T. */
    const char *scaling_env = getenv("RSM_T_SCALING");
    const int scaling_mode = scaling_env && strcmp(scaling_env, "1") == 0;
    const int scaling_T_grid[] = {75, 100, 150, 200, 300, 500, 800, 1200};
    const int scaling_T_len = (int)(sizeof(scaling_T_grid) / sizeof(scaling_T_grid[0]));
    const int num_scen = scaling_mode
        ? scaling_T_len
        : T_len_list_len * rho_list_len * R2_list_len * P;
    Scenario *scen = (Scenario*)malloc(sizeof(Scenario)*num_scen);
    int sidx=0;
    if (scaling_mode) {
        for (int ti=0; ti<scaling_T_len; ti++) {
            scen[sidx++] = (Scenario){scaling_T_grid[ti], 0.8, 0.3,
                                      0.3, 0.7, 0.0};
        }
        printf("T-scaling mode: m=%d rho=0.8 R2=0.3 phis=(0.3,0.7,0.0)\n", N);
    } else {
        for (int ti=0; ti<T_len_list_len; ti++)
        for (int ri=0; ri<rho_list_len; ri++)
        for (int r2=0; r2<R2_list_len; r2++)
        for (int pi=0; pi<P; pi++){
            scen[sidx++] = (Scenario){ T_len_list[ti], rho_list[ri], R2_list[r2], phi_f_arr[pi], phi_eta_arr[pi], phi_u_arr[pi] };
        }
    }

    for (int s=0; s<num_scen; s++){
        const int T_len = scen[s].T_len;
        const double rho = scen[s].rho;
        const double R2_target = scen[s].R2_target;
        const double phi_f = scen[s].phi_f;
        const double phi_eta = scen[s].phi_eta;
        const double phi_u = scen[s].phi_u; /* reported but unused for u_t (iid) */

        const int oos_length = OOS_LEN_FIXED;   // exactly 20 evaluated forecasts
        const int steps = oos_length;
        const int init_window = T_len - oos_length;
        const int start = init_window + 1;
        const int end   = T_len;
        /* One extra terminal observation is needed because Y[t] is predicted
           from X[t-1]. T_len remains the paper's reported sample size. */
        const int sim_len = T_len + 1;

        /* Sigma (and its chol) per scenario */
        double *Sigma = (double*)malloc(sizeof(double)*N*N);
        toeplitz_cov(rho, N, Sigma);
        double *Sigma_chol = (double*)malloc(sizeof(double)*N*N);
        if (chol_copy(Sigma, N, Sigma_chol)!=0){ fprintf(stderr,"Sigma chol failed\n"); return 1; }

        /* Population-R2 scaling based on the stationary covariance
             Var(x_t) = 11'/(1-phi_f^2) + Sigma/(1-phi_eta^2).
           With beta = C*beta_raw and Var(u)=sigma_e^2, choose C so that
             Var(beta'x_t) / (Var(beta'x_t)+sigma_e^2) = R2_target. */
        const double factor_var_scale = 1.0 / (1.0 - phi_f*phi_f);
        const double eta_var_scale = 1.0 / (1.0 - phi_eta*phi_eta);
        double V_signal=0.0;
        for (int i=0;i<N;i++)
            for (int k=0;k<N;k++)
                V_signal += beta_raw[i]
                          * (factor_var_scale + eta_var_scale*Sigma[i*N+k])
                          * beta_raw[k];
        double C = sqrt((R2_target/(1.0-R2_target))
                        * (sigma_e*sigma_e) / V_signal);
        double beta_true[N]; for (int i=0;i<N;i++) beta_true[i] = C*beta_raw[i];
        const double signal_var = C*C*V_signal;
        const double achieved_R2 = signal_var / (signal_var + sigma_e*sigma_e);
        printf("R2 check: target=%.6f analytical=%.6f\n", R2_target, achieved_R2);

        /* z candidates: iterate all from 2..N (inclusive) */
        int z_candidates[N-1]; int zc = 0;
        for (int z = 2; z <= N; z++) z_candidates[zc++] = z;

        /* ===== Baseline ===== */
        double sum_base = 0.0, sum_e_base = 0.0, sum_e2_base = 0.0, cnt_base = 0.0;

        #pragma omp parallel for if(n_rep>1) schedule(dynamic) reduction(+:sum_base,sum_e_base,sum_e2_base,cnt_base)
        for (int repl=0; repl<n_rep; repl++){
            xrng r; norm_state st = {0,0.0};
            unsigned seed = MAIN_SEED ^ (unsigned)(s*1315423911u) ^ (unsigned)(repl*2654435761u);
            xrng_seed(&r, seed);

            double *F_t = (double*)malloc(sizeof(double)*sim_len);
            double *u_t = (double*)malloc(sizeof(double)*sim_len);
            double *Y   = (double*)malloc(sizeof(double)*sim_len);
            double *E_t = (double*)malloc(sizeof(double)*(size_t)sim_len*N);
            double *X   = (double*)malloc(sizeof(double)*(size_t)sim_len*N);
            double *y_tr= (double*)malloc(sizeof(double)*sim_len);
            double *csq = (double*)malloc(sizeof(double)*N);
            double *cxy = (double*)malloc(sizeof(double)*N);
            double *gamma=(double*)malloc(sizeof(double)*N);
            double *f_base=(double*)malloc(sizeof(double)*steps);

            /* DGP */
            gen_AR1_scalar_innovvar(sim_len, phi_f, 1.0, &r, &st, F_t);
            gen_E_AR1_innovcov(sim_len, N, phi_eta, Sigma_chol, &r, &st, E_t);

            for (int t=0;t<sim_len;t++){
                const double Ft = F_t[t];
                double *xrow = X + (size_t)t*N;
                const double *erow = E_t + (size_t)t*N;
                for (int i=0;i<N;i++) xrow[i] = Ft*1.0 + erow[i];
            }

            /* u_t i.i.d. N(0, sigma_e^2) */
            gen_iid_noise(sim_len, sigma_e*sigma_e, &r, &st, u_t);

            Y[0]=NAN;
            for (int t=1;t<sim_len;t++){
                const double *xlag = X + (size_t)(t-1)*N;
                double ssum=0.0; for (int i=0;i<N;i++) ssum += xlag[i]*beta_true[i];
                Y[t] = ssum + u_t[t];
            }

            /* forecast origins */
            for (int k=1; k<=steps; k++){
                int t_now = init_window + k;
                int T_tr = t_now - 1;
                for (int t=0;t<T_tr;t++) y_tr[t] = Y[t+1];
                for (int i=0;i<N;i++){ csq[i]=0.0; cxy[i]=0.0; }
                for (int t=0;t<T_tr;t++){
                    const double *xt = X + (size_t)t*N;
                    const double yt = y_tr[t];
                    for (int i=0;i<N;i++){ double xi=xt[i]; csq[i]+=xi*xi; cxy[i]+=xi*yt; }
                }
                for (int i=0;i<N;i++) gamma[i] = (csq[i]==0.0 ? 0.0 : (cxy[i]/csq[i]));
                const double *x_now = X + (size_t)(t_now-1)*N;
                double acc=0.0; for (int i=0;i<N;i++) acc += gamma[i]*x_now[i];
                f_base[k-1] = acc / (double)N;
            }

            double mse_b=0.0; int cnt=0;
            for (int t=start; t<=end; t++){
                double err = Y[t]-f_base[t-start];
                mse_b += err*err; cnt++;
                sum_e_base  += err;
                sum_e2_base += err*err;
                cnt_base    += 1.0;
            }
            sum_base += sqrt(mse_b / (double)cnt);

            free(F_t); free(u_t); free(Y); free(E_t); free(X);
            free(y_tr); free(csq); free(cxy); free(gamma); free(f_base);
        }

        double mse_base_pooled = (cnt_base>0.0 ? sum_e2_base / cnt_base : NAN);
        double bias2_base = (cnt_base>0.0 ? pow(sum_e_base / cnt_base, 2.0) : NAN);
        double var_base   = fmax(0.0, mse_base_pooled - bias2_base);
        double rmsfe_base = sum_base / (double)n_rep;
        double rmse_pooled_base = sqrt(mse_base_pooled);

        printf("SCEN %d/%d  T=%d rho=%.2f R2=%.2f phis=(%.1f,%.1f,%.1f)\n",
               s+1, num_scen, T_len, rho, R2_target, phi_f, phi_eta, phi_u);
        printf("  Baseline RMSFE(avg): %.6f  |  RMSE(pooled): %.6f  |  bias: %.6f  var: %.6f\n",
               rmsfe_base, rmse_pooled_base, bias2_base, var_base);

        /* ===== RS by z ===== */
        for (int zi=0; zi<zc; zi++){
            const int current_z = z_candidates[zi];

            double sum_rs = 0.0, sum_e_rs = 0.0, sum_e2_rs = 0.0, cnt_rs = 0.0;
            double sum_Q_ols = 0.0, cnt_Q_ols = 0.0;
            double sum_Q_rs_bag  = 0.0, cnt_Q_rs_bag  = 0.0;
            double sum_Q_rs_sub  = 0.0, cnt_Q_rs_sub  = 0.0;
            /* Optional sanity-check losses using I(T,k)=(T-1)/(T-k). */
            double sum_L_ols = 0.0, cnt_L_ols = 0.0;
            double sum_L_rs_bag = 0.0, cnt_L_rs_bag = 0.0;
            double sum_L_rs_sub = 0.0, cnt_L_rs_sub = 0.0;

            #pragma omp parallel for if(n_rep>1) schedule(dynamic) reduction(+:sum_rs,sum_e_rs,sum_e2_rs,cnt_rs,sum_Q_ols,cnt_Q_ols,sum_Q_rs_bag,cnt_Q_rs_bag,sum_Q_rs_sub,cnt_Q_rs_sub,sum_L_ols,cnt_L_ols,sum_L_rs_bag,cnt_L_rs_bag,sum_L_rs_sub,cnt_L_rs_sub)
            for (int repl=0; repl<n_rep; repl++){
                xrng r; norm_state st = {0,0.0};
                unsigned seed = MAIN_SEED ^ (unsigned)(s*1315423911u) ^ (unsigned)(repl*2654435761u);
                xrng_seed(&r, seed);

                double *F_t = (double*)malloc(sizeof(double)*sim_len);
                double *u_t = (double*)malloc(sizeof(double)*sim_len);
                double *Y   = (double*)malloc(sizeof(double)*sim_len);
                double *E_t = (double*)malloc(sizeof(double)*(size_t)sim_len*N);
                double *X   = (double*)malloc(sizeof(double)*(size_t)sim_len*N);
                double *y_tr= (double*)malloc(sizeof(double)*sim_len);
                double *csq = (double*)malloc(sizeof(double)*N);
                double *cxy = (double*)malloc(sizeof(double)*N);
                double *gamma=(double*)malloc(sizeof(double)*N);
                double *S_lower = (double*)malloc(sizeof(double)*N*N);
                double *b = (double*)malloc(sizeof(double)*N);
                double *s1= (double*)malloc(sizeof(double)*N);
                int *pool_idx = (int*)malloc(sizeof(int)*N);
                int *sub_idx  = (int*)malloc(sizeof(int)*N);
                double *S_sub = (double*)malloc(sizeof(double)*N*N);
                double *f_rs  = (double*)malloc(sizeof(double)*steps);
                double *Sigma_u_lower = (double*)malloc(sizeof(double)*N*N);

                /* DGP */
                gen_AR1_scalar_innovvar(sim_len, phi_f, 1.0, &r, &st, F_t);
                gen_E_AR1_innovcov(sim_len, N, phi_eta, Sigma_chol, &r, &st, E_t);
                for (int t=0;t<sim_len;t++){
                    const double Ft = F_t[t];
                    double *xrow = X + (size_t)t*N;
                    const double *erow = E_t + (size_t)t*N;
                    for (int i=0;i<N;i++) xrow[i] = Ft*1.0 + erow[i];
                }
                /* u_t iid */
                gen_iid_noise(sim_len, sigma_e*sigma_e, &r, &st, u_t);

                Y[0]=NAN;
                for (int t=1;t<sim_len;t++){
                    const double *xlag = X + (size_t)(t-1)*N;
                    double ssum=0.0; for (int i=0;i<N;i++) ssum += xlag[i]*beta_true[i];
                    Y[t] = ssum + u_t[t];
                }

                /* forecast origins */
                for (int k=1; k<=steps; k++){
                    int t_now = init_window + k;
                    int T_tr = t_now - 1;      /* per-origin training size */
                    for (int t=0;t<T_tr;t++) y_tr[t] = Y[t+1];

                    /* gamma */
                    for (int i=0;i<N;i++){ csq[i]=0.0; cxy[i]=0.0; }
                    for (int t=0;t<T_tr;t++){
                        const double *xt = X + (size_t)t*N;
                        const double yt = y_tr[t];
                        for (int i=0;i<N;i++){ double xi=xt[i]; csq[i]+=xi*xi; cxy[i]+=xi*yt; }
                    }
                    for (int i=0;i<N;i++) gamma[i] = (csq[i]==0.0 ? 0.0 : (cxy[i]/csq[i]));

                    /* grams (kept for compatibility) */
                    double Ty;
                    build_grams(X, y_tr, T_tr, N, S_lower, b, s1, &Ty);

                    const double *x_now = X + (size_t)(t_now-1)*N;

                    /* --- Build Σ_u on training window --- */
                    build_sigma_u_lower(X, y_tr, gamma, T_tr, N, Sigma_u_lower);

                    /* --- Compute full-information Q_ols & L_ols at this origin --- */
                    {
                        double w_star_tmp[128];
                        double QN = NAN;
                        if (sum_to_one_weights_from_cov_lower(Sigma_u_lower, N, w_star_tmp, &QN)==0){
                            double infl_ols = (T_tr > N)
                                ? ((double)T_tr - 1.0) / ((double)T_tr - (double)N)
                                : NAN;
                            sum_Q_ols += QN;          cnt_Q_ols += 1.0;
                            if (isfinite(infl_ols)) {
                                sum_L_ols += QN * infl_ols; cnt_L_ols += 1.0;
                            }
                        }
                    }

                    /* base forecasts at now */
                    double f_now[128]; for (int i=0;i<N;i++) f_now[i] = gamma[i]*x_now[i];

                    double pred_k = 0.0;

                    if (current_z == N){
                        /* Use full-pool constrained-optimal weights when z=N.
                           Populate the RSM Q/L columns too, because at this
                           endpoint RSM, subset-average, and full-pool coincide. */
                        double wN[128]; double QNdummy = NAN;
                        if (sum_to_one_weights_from_cov_lower(Sigma_u_lower, N, wN, &QNdummy)==0){
                            double acc=0.0; for (int i=0;i<N;i++) acc += wN[i]*f_now[i];
                            pred_k = acc;
                            sum_Q_rs_bag += QNdummy; cnt_Q_rs_bag += 1.0;
                            sum_Q_rs_sub += QNdummy; cnt_Q_rs_sub += 1.0;
                            if (T_tr > N) {
                                double infl_full = ((double)T_tr - 1.0) /
                                                   ((double)T_tr - (double)N);
                                sum_L_rs_bag += QNdummy * infl_full; cnt_L_rs_bag += 1.0;
                                sum_L_rs_sub += QNdummy * infl_full; cnt_L_rs_sub += 1.0;
                            }
                        } else pred_k = 0.0;
                    } else {
                        /* RS: average predictions; build bagged weight; average per-subset Q */
                        double acc_rs = 0.0;
                        double w_bar[128]; for (int i=0;i<N;i++) w_bar[i]=0.0;
                        double acc_Q_sub = 0.0;
                        int succ = 0;

                        for (int rr=0; rr<K_RS; rr++){
                            sample_k_no_replace(&r, N, current_z, pool_idx, sub_idx);
                            gather_lower(Sigma_u_lower, N, sub_idx, current_z, S_sub);

                            double wS[128]; double QS = NAN;
                            if (sum_to_one_weights_from_cov_lower(S_sub, current_z, wS, &QS)==0){
                                double pr = 0.0;
                                for (int j=0;j<current_z;j++){
                                    int idx = sub_idx[j];
                                    pr += wS[j] * f_now[idx];
                                    w_bar[idx] += wS[j];
                                }
                                acc_rs   += pr;
                                acc_Q_sub += QS;
                                succ++;
                            }
                        }

                        if (succ > 0){
                            pred_k = acc_rs / (double)succ;
                            for (int i=0;i<N;i++) w_bar[i] /= (double)succ;
                            double sumw=0.0; for (int i=0;i<N;i++) sumw += w_bar[i];
                            if (fabs(1.0 - sumw) > 1e-10){
                                for (int i=0;i<N;i++) w_bar[i] /= sumw;   /* enforce 1'w=1 */
                            }
                            double Q_bag = quadform_lower(Sigma_u_lower, N, w_bar);
                            double Q_subavg = acc_Q_sub / (double)succ;

                            double infl_rs = (T_tr > current_z)
                                ? ((double)T_tr - 1.0) /
                                  ((double)T_tr - (double)current_z)
                                : NAN;

                            sum_Q_rs_bag += Q_bag;           cnt_Q_rs_bag += 1.0;
                            sum_Q_rs_sub += Q_subavg;        cnt_Q_rs_sub += 1.0;
                            if (isfinite(infl_rs)) {
                                sum_L_rs_bag += Q_bag * infl_rs; cnt_L_rs_bag += 1.0;
                                sum_L_rs_sub += Q_subavg * infl_rs; cnt_L_rs_sub += 1.0;
                            }
                        } else {
                            pred_k = 0.0; /* extremely rare */
                        }
                    }

                    f_rs[k-1] = pred_k;
                } /* origins */

                /* errors for this z */
                double mse=0.0; int c=0;
                for (int t=start; t<=end; t++){
                    double err = Y[t]-f_rs[t-start];
                    mse += err*err; c++;
                    sum_e_rs  += err;
                    sum_e2_rs += err*err;
                    cnt_rs    += 1.0;
                }
                sum_rs += sqrt(mse / (double)c);

                free(F_t); free(u_t); free(Y); free(E_t); free(X);
                free(y_tr); free(csq); free(cxy); free(gamma);
                free(S_lower); free(b); free(s1);
                free(pool_idx); free(sub_idx); free(S_sub);
                free(f_rs); free(Sigma_u_lower);
            } /* repl loop */

            double avg_rs = sum_rs / (double)n_rep;
            double rel = avg_rs / rmsfe_base;

            double mse_rs_pooled = (cnt_rs>0.0 ? sum_e2_rs / cnt_rs : NAN);
            double bias2_rs = (cnt_rs>0.0 ? pow(sum_e_rs / cnt_rs, 2.0) : NAN);
            double var_rs   = fmax(0.0, mse_rs_pooled - bias2_rs);
            double rmse_pooled_rs = sqrt(mse_rs_pooled);

            double avg_Q_ols     = (cnt_Q_ols    >0.0 ? (sum_Q_ols    / cnt_Q_ols    ) : NAN);
            double avg_Q_rs_bag  = (cnt_Q_rs_bag >0.0 ? (sum_Q_rs_bag / cnt_Q_rs_bag ) : NAN);
            double avg_Q_rs_sub  = (cnt_Q_rs_sub >0.0 ? (sum_Q_rs_sub / cnt_Q_rs_sub ) : NAN);

            double avg_L_ols     = (cnt_L_ols    >0.0 ? (sum_L_ols    / cnt_L_ols    ) : NAN);
            double avg_L_rs_bag  = (cnt_L_rs_bag >0.0 ? (sum_L_rs_bag / cnt_L_rs_bag ) : NAN);
            double avg_L_rs_sub  = (cnt_L_rs_sub >0.0 ? (sum_L_rs_sub / cnt_L_rs_sub ) : NAN);

            printf("  z=%2d: RS=%.6f  rel=%.6f%s  |  RMSE(pooled): %.6f  |  bias: %.6f  var: %.6f  |  Q_ols: %.6f  Q_rs(bag): %.6f  Q_rs(subavg): %.6f  |  L_ols: %.6f  L_rs(bag): %.6f  L_rs(subavg): %.6f\n",
                   current_z, avg_rs, rel, (current_z==N ? " [sum-to-one full]" : ""),
                   rmse_pooled_rs, bias2_rs, var_rs,
                   avg_Q_ols, avg_Q_rs_bag, avg_Q_rs_sub,
                   avg_L_ols, avg_L_rs_bag, avg_L_rs_sub);

            fprintf(fp,
                    "%d,%.6f,%.6f,%.1f,%.1f,%.1f,%d,"
                    "%.6f,%.6f,%.6f,%.6f,%.6f,"
                    "%.6f,%.6f,%.6f,%.6f,"
                    "%.6f,%.6f,%.6f,"      /* Q_ols, Q_rs_bag, Q_rs_sub */
                    "%.6f,%.6f,%.6f\n",    /* L_ols, L_rs_bag, L_rs_sub */
                    T_len, rho, R2_target, phi_f, phi_eta, phi_u,
                    current_z, rmsfe_base, avg_rs, rel,
                    rmse_pooled_base, rmse_pooled_rs,
                    bias2_base, var_base, bias2_rs, var_rs,
                    avg_Q_ols, avg_Q_rs_bag, avg_Q_rs_sub,
                    avg_L_ols, avg_L_rs_bag, avg_L_rs_sub);


            fflush(fp);
        } /* z loop */

        free(Sigma); free(Sigma_chol);
    }

    free(scen);
    fclose(fp);
    printf("\nSaved results to rs_results.csv\n");
    return 0;
}

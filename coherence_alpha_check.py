import numpy as np
from itertools import combinations
def alpha_exact(tau):
    m=len(tau); tb=tau.mean(); ks=range(1,m); R=[]
    for k in ks:
        num=np.zeros(m); cnt=np.zeros(m)
        for S in combinations(range(m),k):
            A=tau[list(S)].sum()
            for i in S: num[i]+=1/A; cnt[i]+=1
        g=k*num/cnt; R.append((1/m**2)*np.sum(tau*(g-1/tb)**2))
    z=np.array([1/k-1/m for k in ks])
    return np.polyfit(np.log(z),np.log(R),1)[0]
np.random.seed(1)
rows=[]
profiles={
 "homog-ish":1+0.05*np.random.randn(9),
 "mild":1+0.2*np.random.randn(9),
 "lognormal0.5":np.exp(0.5*np.random.randn(9)),
 "two-group":np.array([5.,5,1,1,1,1,1,1,1]),
 "lognormal1.0":np.exp(1.0*np.random.randn(9)),
 "dominant":np.array([12.,1,1,1,1,1,1,1,1]),
 "very-dom":np.array([40.,1,1,1,1,1,1,1,1]),
}
print(f"{'profile':14s} {'coherence kappa':>15s} {'alpha':>8s}")
for name,t in profiles.items():
    t=np.abs(t)
    kappa=t.max()/t.mean()
    print(f"{name:14s} {kappa:15.2f} {alpha_exact(t):8.3f}")

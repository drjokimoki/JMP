# Theory review: current `Job_Market_paper_Kozyrev-2.pdf`

Reviewed file: 61 pages, generated July 11, 2026 at 15:38 CEST. Scope: Sections 3--4 and Appendices A--G only.

## Verdict

The exact population algebra is largely valid. The current PDF is improved relative to the previous version: it now says that Equation (4) uses fixed weights, identifies the inflation factor as an approximation for overlapping subset regressions, removes the duplicated average-subset calculation from Appendix B, and weakens the earlier monotonic-interpolation wording.

The theory still needs several corrections. The most important is the feasible domain in Proposition 8. The rate result also needs to be identified explicitly as asymptotic guidance for a sequence of finite-sample criteria, rather than an exact finite-sample property of the implemented estimator.

## Substantive corrections

### 1. Proposition 8 uses the wrong feasible set

The criterion

\[
L_k^{\mathrm{bag}}=(A+R_k^{\mathrm{bag}})\frac{T-1}{T-k}
\]

requires `k<T`. The paper instead defines the minimizer over `1 <= k <= m`. If `m >= T`, this domain contains `k=T`, where the criterion is undefined, and `k>T`, where the inflation factor becomes negative. Replace the domain everywhere in Proposition 8 and Appendix G by

\[
1\le k\le K_T,\qquad K_T=\min\{m,T-1\}.
\]

This is a genuine logical error, not merely a presentation issue.

### 2. Separate finite-sample criteria from asymptotic rate statements

The criterion is evaluated for finite `T,m,k`, but the notations `O(k^{-3})`, `Omega(T^{1/3})`, and `O(T^{1/2})` refer to sequences. There is no contradiction, but the sequence must be stated. For Proposition 8, write `m=m_T`, let `T -> infinity`, and impose `T^{1/2}/m_T -> 0`. The constants in Assumption 4 must be uniform along that sequence, not merely “independent of T and k.”

The result should be titled “Rate window under the stylized finite-sample criterion.” It is not an exact rate theorem for the feasible finite-`G` estimator because overlapping regressions have correlated coefficient-estimation errors.

### 3. Proposition 5 is an approximation unless an asymptotic sequence is supplied

The sentence “if the precisions are bounded away as mentioned in Assumption 3” is incorrect: Assumption 3 does not impose a uniform lower and upper bound on the sequence of precision profiles. For a finite-sample presentation, use

\[
E_S[1/A_S]\approx \frac1{k\bar\tau}
+\frac{m-k}{k^2(m-1)}\frac{v_\tau}{\bar\tau^3}
\]

and describe it as a second-order finite-population Taylor approximation. If the `O(k^{-3})` remainder is retained, explicitly introduce a sequence with `m -> infinity`, `k -> infinity`, and `0<underline tau <= tau_{i,m} <= overline tau < infinity` uniformly. For the later rate argument, the stronger interior regime `k/m -> 0` is needed.

### 4. The joint moment display is dimensionally wrong and overstates distributional assumptions

The mean of `f_t` is `mu_t iota_m`, not the scalar `mu_t`. Moreover, the symbol `sim` suggests a fully specified joint distribution, although Assumptions 1--2 specify only moments. Replace the display by explicit conditional mean and covariance equations. State whether expectations are conditional on the forecast information set; otherwise “unbiased conditional forecasts” is stronger than `E(u_t)=0`.

### 5. Equation (4) needs the full fixed-weight qualification in the text

The new footnote is helpful but too short. Explain that population weights are fixed and that estimated weights are evaluated on a new observation independent of the training sample. Conditional on the training sample, estimated weights are fixed; unconditional risk then averages over training-sample uncertainty. Without sample separation, additional covariance terms can occur.

### 6. Proposition 2 still overclaims after the endpoint result

The endpoint result is valid: `k=1` gives equal weights and `k=m` gives population-optimal weights. But the following paragraph still says that RSM “detects” homogeneity, “selects” `k=m`, and will “presumably get results not worse than any benchmark.” None follows from Proposition 2. Data-dependent selection of `k` has its own tuning cost. Replace this paragraph with the endpoint-only interpretation in the accompanying LaTeX file.

### 7. Proposition 3's interpretation remains too strong and grammatically broken

The sentence “cannot obtain a weight exceeding a equal weight having an above 1/k conditional mean...” is missing “without” and is difficult to parse. More importantly, the result concerns a conditional mean weight; it does not say that the forecast is useful or accurate in each subset. Use the replacement text supplied below.

### 8. Distinguish population bagging from estimated RSM

Proposition 1 is correct for fixed population subset-optimal weights. `Q_k^{bag,G}` is therefore a population risk averaged over random subset draws, not the actual risk “a decision-maker gets” after estimating the subset coefficients. State this immediately after Definition 2. The exact `1/G` bridge does not include training-sample coefficient error.

## Results that are valid as written or after minor notation fixes

- Proposition 1: exact for independent subset draws and fixed population subset weights.
- Proposition 2: both endpoint identities.
- Proposition 3: the algebraic equivalence under uniform fixed-size sampling.
- Proposition 4: inverse-precision subset weights and `Q(S)=q+1/A_S`.
- Section 4.4: optimal/equal-weight risks and the `B,H` comparison.
- Proposition 6: exact oracle-distance representation.
- Proposition 7: exact Jensen ceiling and projection floor, provided the heterogeneous case makes the denominator positive.
- Proposition 8: its proof is valid conditional on the corrected feasible domain, a triangular asymptotic sequence, uniform Assumption 4, and the stipulated loss criterion.

## Duplication and organization

- The earlier Appendix B/C duplication has been fixed.
- Appendix D still introduces “Proposition 9,” although the same square-root result already appears in Section 4.5. Number the result only once. Appendix D should be “Proof of the average-subset square-root approximation.”
- Section 4.6 states Proposition 6 and Appendix E appropriately proves it; keep this pattern consistently for all propositions.
- Section 3 and Section 4.1 still repeat parts of the random-subset procedure. Keep Section 3 algorithmic and Section 4 devoted to population risks.

## Presentation corrections

- Use “random-subset method” consistently, including capitalization in the Section 4.1 heading.
- In Section 3, call `R_i` a selection matrix, not a permutation matrix; it is rectangular.
- Harmonize `N,m,n,T`. In Section 3 the target is incorrectly called an `n`-dimensional vector even though it is scalar at each date.
- Use one notation for precision variance: `v_tau`, `V_tau`, or `V_delta`, and define it once.
- Replace `Sigma_{e.m}` by `Sigma_{e,m}`.
- Use `K_T=min{m,T-1}` in every practical clipping rule.
- Add `\hypersetup{hidelinks}` to remove visible red boxes around references in the rendered PDF.

Paste-ready corrections are in `theory_corrections_current_pdf.tex`.

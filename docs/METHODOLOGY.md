# Methodology, conventions and design decisions

Notation follows the specification: `theta` = mean reversion speed, `mu` = long-term mean, `sigma` = volatility,
`r0` = initial rate (Java: `meanReversionSpeed`, `longTermMean`, `volatility`, `initialRate`). `C_i` are bond
cashflows (`bondCashflow`), `C_j` / `\hat C_j` continuation values (`continuationValue`). Times are in years.

## 1. Model and contract

**CIR under Q.** `dr = theta (mu - r) dt + sigma sqrt(r) dW`, all four parameters `> 0` in both Feller modes.
`STRICT` (default) rejects `2 theta mu < sigma^2` when the parameters are created (`FellerConditionViolationException`);
`ALLOW_VIOLATION` must be chosen explicitly. The check is done **exactly in decimal arithmetic** on the input
strings (e.g. `theta = 0.5, mu = 0.04, sigma = 0.2` lies exactly on the boundary and is accepted, whereas binary64
would compute `0.2 * 0.2 > 0.04`). The Feller condition concerns the continuous process (origin not attained); it
does not prevent negative values of the Euler approximation (section 4).

**Inputs.** All business inputs are decimal strings. BigDecimal values are parsed from them; double values are
`Double.parseDouble` of the same strings (correctly rounded). Both variants are thus derived from one configuration.

**Underlying.** Positive cashflows `C_1..C_n` at `0 < t_1 < ... < t_n`; `CashflowSchedule.couponBond` creates the
regular bond `C_1 = ... = C_{n-1} = cN`, `C_n = (1+c)N` (frequency `f`: `cN/f` at `k/f`, exact decimal division
required). A single cashflow is allowed (European / zero-bond cases).

**Option.** Bermudan call with money strike `K > 0` (same unit as the cashflows, no implicit percentage or notional
scaling) and exercise dates `0 < e_1 < ... < e_m < t_n`. A single exercise date is the European option on the same
underlying.

**Binding cashflow convention.** At `e_j` the holder buys the cashflows with `t_i > e_j` strictly; a cashflow paid
exactly at `e_j` is not part of the remaining bond. `B_{e_j} = sum_{t_i > e_j} C_i P(e_j, t_i)`,
`E_j = (B_{e_j} - K)^+` (immediate exercise value, distinct from the continuation value).

**Zero bonds.** Analytic CIR formula, never simulated:

```
P(t,T) = A(tau) exp(-B(tau) r_t), tau = T - t, h = sqrt(theta^2 + 2 sigma^2)
B(tau) = 2 (e^{h tau} - 1) / (2h + (theta + h)(e^{h tau} - 1))
A(tau) = [2h e^{(theta+h) tau/2} / (2h + (theta + h)(e^{h tau} - 1))]^{2 theta mu / sigma^2}
```

(Cox, Ingersoll, Ross 1985; Brigo & Mercurio, *Interest Rate Models - Theory and Practice*, 2nd ed. 2006, sect. 3.2.3.)
It is evaluated as `exp(ln A - B r)` with `expm1` for `e^{h tau} - 1` to avoid overflow of the power for large
exponents `c = 2 theta mu / sigma^2` (small sigma) and cancellation for short maturities. `ln A` and `B` depend only on
`t_i - e_j` and are precomputed per (exercise date, cashflow). Two algebraically equivalent forms of `ln A` exist
(`zerobond.form`):

* `LOG_DIFFERENCE` (textbook, default): `c (ln 2h + (theta+h) tau/2 - ln denom)`. The bracket is a small
  difference of O(1) terms; in binary64 this costs up to ~2e-15 relative accuracy of P for tau of a few years
  (measured against mpmath) and ~4e-12 for `c = 128, tau = 30`. **This is the dominant source of
  double-vs-BigDecimal differences in the exercise values** (see RESULTS).
* `LOG1P_RATIO`: `c log1p((2h expm1(g) - (theta+h) em1) / denom)` with `g = (theta+h) tau/2`; better for short and
  medium maturities, for very long maturities the numerator cancels instead.

Both forms are implemented in both arithmetics with the same operation order.

## 2. Backward induction

On every training path the induction keeps **only the selected later exercise index `tau^p`** (or -1), not a tree of
intermediate values; cashflows and discount factors are recomputed from the stored exercise values and cumulative
integrals.

```
last date m:  exercise iff E_m > 0                                     (V_m = E_m)
j = m-1..1:   Y_j^p = exp(-(I_{tau^p} - I_{e_j})) E_{tau^p}             (0 if tau^p = -1: no later exercise)
              \hat C_j^p = estimate from {Y_j^q : q in estimation set}
              exercise iff E_j^p > 0 and E_j^p > \hat C_j^p           (tie E = C: continue; never exercise E = 0)
              exercise: tau^p := j (replaces the later cashflow)
price:        V_0 = mean_p exp(-I_{tau^p}) E_{tau^p}
```

The decisions stored per date are *hypothetical* decisions of the induction; the realised exercise of a path is its
first positive decision. Training (policy estimation) and valuation (applying the frozen policy forward, in
chronological order, to independent paths) are separated (section 7).

## 3. Continuation estimation

### Quantile clusters (primary)

State variable: effective short rate `x^+_{e_j}` (CIR is one-dimensional and Markov).

* Sort the estimation set by `(state, pathIndex)`; the path index is the deterministic tie-breaker.
* Cluster `k` = sorted positions `[floor(k n/K), floor((k+1) n/K))`; sizes differ by at most one.
* `\hat C_{j,k} = (1/|G_k|) sum_{p in G_k} Y_j^p`, summed in sorted order; piecewise constant in the state.
* Minimum size: `K_eff = max(1, min(K, floor(n / minPathsPerCluster)))` (default min 20); fewer eligible paths than
  the minimum -> one cluster. Recommendation: at least 50-100 paths per cluster for a standard error of the cluster
  mean below ~10% of the cluster dispersion; `K` 10-50 for 10^3..10^5 paths. Too few clusters -> approximation bias
  (coarse step function), too many -> variance of the local means (and noisy boundary).
* Ties: identical states can be split between two neighbouring clusters (by path index). They occur in practice
  only for `x^+ = 0` after truncation (Feller-violating runs).
* Boundaries for new states: `u_k` = state of the last training path of cluster `k`; cluster `k` covers
  `(u_{k-1}, u_k]`, states below/above the training range belong to the first/last cluster.

**Bond value as cluster variable?** `B_{e_j}(r)` is strictly decreasing in `r`, so sorting by `B` gives the same
order (reversed) and hence the *same* quantile clusters (up to ties). Clustering on `B` is therefore equivalent for
one factor; it would matter only for multi-factor models or for non-quantile (equal-width) clusters. For LSM the
choice is not equivalent (different basis), see below.

### Longstaff-Schwartz (second method)

Same targets `Y_j^p`; basis `1, x, x^2` with `x = r` (default) or `x = B` (`lsm.basisVariable`). Solved by
Householder QR (`Q^T y`, back substitution, no normal equations, no explicit inverse) in both arithmetics with the
same operation order; rank-deficient columns get coefficient 0; fewer than 3 paths -> constant (mean) estimate.
No scaling of the basis is applied on purpose: `1, B, B^2` with `B ~ 100` is badly conditioned (cond ~1e7..1e8),
which makes the precision effect of the linear algebra visible (experiment `sens-lsmbasis`). For production one would
centre/scale the variable.

### Estimation set (decision, documented by tests)

Default `ALL`: all training paths at `e_j` (also out-of-the-money ones, whose targets carry information about the
continuation value near the money). Exercise is still restricted to `E_j > 0`. `IN_THE_MONEY` (classic LSM) is
available. In both cases the selection depends only on `E_j` (state at `e_j`), never on later-computed earlier
decisions. In Mode A, the estimation set is the one of the double variant for both arithmetics; differences in the
selection (Mode B, ITM) appear as cluster switches. Paths outside the estimation set are never exercised.

## 4. Simulation, Euler and negative values

Full-truncation Euler on a uniform grid that contains every exercise date exactly (checked in decimal):

```
x_{k+1} = x_k + theta (mu - x_k^+) dt + sigma sqrt(x_k^+) dW_k,   x^+ = max(x, 0),   x_0 = r0
```

The internal state `x_k` may become negative and is **not** projected (projection after every step would be a
different scheme). The effective rate `x_k^+` is used for drift, diffusion, the discount integral, bond values,
clustering and regression. Negative internal states are counted per path in both arithmetics (`Euler` line of the
report). Negative Euler values are a discretisation artefact; they can occur even when the Feller condition holds
(probability per step ~ P(Z < -x/(sigma sqrt(x dt))), negligible for the standard case: 0 occurrences in 1000 paths,
but frequent for sigma = 0.3-0.4). Milstein and exact (noncentral chi-square) transitions are later options;
exact transitions at grid points do not make the integrated rate path exact.

The order of floating-point operations is mirrored line by line in the BigDecimal generator:
`x + (theta (mu - x^+)) dt + (sigma sqrt(x^+)) dW`.

## 5. Discounting

Left rectangle rule on the simulation grid: `I_{k+1} = I_k + x_k^+ dt`; the cumulative integral is stored at every
exercise date and `D(e_j, e_l) = exp(-(I_l - I_j))`. Same rule and same effective state in both arithmetics. The
rectangle rule has O(dt) bias; the trapezoidal rule is a later control variant. The discretisation bias of the
(Euler + rectangle) pair is measured against analytic and FD references (RESULTS).

## 6. Random numbers and common random numbers

* Generator: per path `p` a `SplittableRandom` (SplitMix64; algorithm stable since Java 8) seeded with
  `mix64(seed + 0x9E3779B97F4A7C15 (p+1))`; uniforms `u = (k + 0.5) 2^-52` from the top 52 bits (exact, strictly in
  (0,1)); normals by Wichura's AS 241 (relative accuracy ~1e-16, validated against mpmath).
* Path `p` always gets the same normals, independently of the number of paths, call order or threads; nothing is
  stored (memory-friendly blocks = one path at a time). Both arithmetics call the same source -> identical normals
  (test `bothArithmeticsConsumeTheSameNormals`).
* BigDecimal receives `new BigDecimal(z)` (exact binary value). The MathContext is applied first in
  `sqrt(dt) * Z`. `BigDecimal.valueOf(z)` (shortest decimal string) would be a different convention and is not used.
  Random inputs are therefore limited to binary64 precision in both variants - BigDecimal cannot recover information
  that the normals never had.
* Brownian coupling: `simulation.normalsPerStep = m` uses `dW = sqrt(dt/m)(Z_1 + ... + Z_m)`; grids with the same
  `stepsPerYear * m` share the same fine Brownian path (e.g. 4, 12, 36, 84, 252 steps on a 252/year base).
* Training and valuation sets use different seeds; within each set both arithmetics share the normals.

## 7. Comparison modes, out-of-sample protocol

* **Mode A (shared clusters).** The BigDecimal induction uses the *training structure built from the double
  states*: identical eligible path indices and identical cluster membership; only the arithmetic of targets,
  discounting, sums and means differs. For valuation paths, Mode A assigns clusters with the double boundaries applied
  to the double valuation states, so assignment is identical as well.
* **Mode B (separate clusters).** Each arithmetic sorts its own states and builds its own clusters and boundaries.
  Small state differences may change sort order, membership, means, decisions.
* The double result is the same in both modes. Shared training clusters do not imply identical policies (decisions
  use each arithmetic's own `E` and `\hat C`). The **frozen-policy control** applies the double policy's realised
  exercise dates to the BigDecimal valuation paths: it isolates pure cashflow arithmetic from decision effects.
* **Out-of-sample:** policies (cluster bounds + means, or coefficients) are trained on training paths and applied in
  chronological order to independent valuation paths. The out-of-sample value of an admissible policy is, in
  expectation, a lower bound of the optimal value *of the simulated (discretised) model*; individual estimates and
  the continuous model may differ.

## 8. Precision strategy

* Working precision: `MathContext` (default 34 digits HALF_EVEN = DECIMAL128), applied to every add/subtract/
  multiply/divide of the valuation and to the results of elementary functions. Configurable precision and rounding
  mode; reference runs at 64 and 100 digits.
* Elementary functions: `ch.obermuhlner:big-math` 2.3.2 (`exp`, `log`, `sqrt`) evaluated with `precision + 10` guard
  digits, then rounded to the working precision; `expm1` adds the number of cancelling leading digits. No
  `Math.exp(x.doubleValue())` anywhere in the BigDecimal valuation.
* No premature rounding to fixed decimal places; output rounding happens only in the console report (CSV: full
  decimal strings).
* Differences `decimal - double` are computed in BigDecimal using the exact value `new BigDecimal(double)` and only
  the (tiny) difference is converted to double for aggregated statistics (max/mean/RMSE/ULP).
* Explicit exceptions / boundaries of the "pure BigDecimal" variant:
  * normals come from a binary64 source (common random numbers);
  * sorting compares BigDecimal values exactly (`compareTo`), ties by path index;
  * the analytic European/coupon-bond references and the FD reference are **double-based** and labelled as such
    (the noncentral chi-square CDF uses Commons Math's regularized gamma function; validated range below);
  * diagnostics (boundary brackets, ULP statistics) convert BigDecimal values to double for reporting only.
* A higher-precision run is a numerical reference **of the same algorithm**, not a proof that the valuation method is
  correct.

### Where a full BigDecimal implementation is problematic

| area | issue | handling |
|---|---|---|
| sqrt | no `BigDecimal.sqrt` accuracy control before Java 9, cost | big-math with guard digits |
| exp/log | not in the JDK; cost grows ~quadratically with digits | big-math with guard digits, precomputed `ln A`, `B` |
| normals | inverse CDF in decimal would be expensive and changes the random input | shared binary64 normals, exact conversion |
| sorting/cluster bounds | exact comparison; different values than double -> different order (Mode B) | Mode A/B separation, switch diagnostics |
| linear algebra (LSM) | square roots and divisions at working precision; conditioning | Householder QR in both arithmetics |
| noncentral chi-square | high-precision CDF is a project of its own | double reference, labelled, validated range |

## 9. Metrics (section 12 of the specification)

| spec | metric | where |
|---|---|---|
| 12.1 | price, SE, 95% CI, runtime per variant; `Delta V`, rel., `|Delta V|/SE` | console, `summary.csv`, `comparisons.csv` |
| 12.2 | short-rate differences at exercise dates (all paths) and along the grid (traced paths): max/mean abs, max rel, ULP | `state_differences.csv`, `trajectory_differences.csv` |
| 12.3 | discount factor `exp(-I)` differences | `state_differences.csv` |
| 12.4 | per date and cluster: `|C_dec - C_dbl|`, rel., sizes, state range, target s.d., `|dC| / (s/sqrt(n))` | `clusters.csv` |
| 12.5 | Mode B cluster switches; rank distance of switched paths to the nearest quantile boundary; `|dr|` | `cluster_switches.csv` |
| 12.6 | per date: agree/differ, double-ex/dec-cont, dec-ex/double-cont (in-sample and out-of-sample) | `exercise_dates.csv` |
| 12.7 | `|E - \hat C|` of paths with differing decisions: median/mean/max (each variant) | `exercise_dates.csv` |
| 12.8 | paths with different later exercise index before the decision at `e_j`; `max |dY_j|`; later decision differences | `propagation.csv` |
| 12.9 | boundary bracket `[max exercised rate, min continued in-the-money rate]` and overlap count | `exercise_dates.csv` |
| 12.10 | `V_LSM - V_Cluster` vs `V_dec - V_dbl`, realised exercise rates, boundaries, differing realised dates | console |
| 12.11 | paired `d_p = Y0_dec - Y0_dbl`: mean, `SE_d`, 95% CI, max, RMSE, paths with `d != 0`, paths with different exercise date | `comparisons.csv` |
| 12.12 | out-of-sample valuation with frozen policies; realised vs hypothetical decisions | all of the above with `sample=out-of-sample` |

Note: `Delta V = V_dec - V_dbl` uses the double mean (computed in double), whereas the paired mean difference is the
exact mean of `d_p`. Their difference is the rounding error of the double aggregation itself.

The paired CI is conditional on the trained policies; training uncertainty is assessed by repeating with other
training seeds (`sens-seeds`). Relative differences use only references above a floor; references at or below the
floor are counted separately (`near_zero_refs`).

**Error categories.** Numerical (arithmetic) effects are `Delta V`, `d_p`, state/discount/cluster differences;
statistical effects are SE/CI; time-discretisation effects are measured by Euler step series against analytic/FD
references; continuation-approximation effects by in-sample vs out-of-sample, cluster count, cluster vs LSM and
MC vs FD. A distance to a reference generally contains several of these at once.

## 10. Analytic and deterministic references

* **Zero bonds** (double and BigDecimal): validated against mpmath (60 digits, decimal inputs) including very short
  maturity (1e-6), small volatility (`c = 40000`), long maturity with `c = 128`, near-Feller parameters. Error model
  used in the tests: relative error `<= 2e-15 (1 + c max(1, tau))` in double and `1e-32 (1 + c max(1, tau))` at 34
  digits.
* **Noncentral chi-square CDF** (double): Poisson mixture summed from the mode, weights normalised by the mass used,
  Commons Math regularized gamma; validated against a 40-digit mpmath series: relative error < 1e-12 for
  noncentrality <= 1000 (all option cases used here have lambda ~ 10-200), ~1e-11 at lambda ~ 1e4 (Commons Math's
  gamma function at a ~ 6000 limits it; SciPy is 1e-14 there).
* **Zero-bond call** closed form and **Jamshidian decomposition** for the coupon-bond call: critical rate by
  bisection on `[0, r_hi]`; if `K >= B(r = 0) = sum c_i A_i` there is no root and the call is worth 0 (rates are
  non-negative). Validated against an independent SciPy implementation.
* **Finite differences** for the Bermudan (and European) option: CIR PDE in time to go, Crank-Nicolson with
  Rannacher start-up after each exercise date, three-point differences on a non-uniform grid, upwind drift only
  where a central difference would produce a negative off-diagonal (near r = 0), PDE itself at r = 0 (diffusion
  vanishes, inflow drift), `V_rr = 0` at `r_max`, `r0` on a node. Default: sinh-stretched grid
  `r = r0 + w sinh(xi)` with `w = max(mu, r0)/2`, and **kink smoothing**: at the final payoff and at each exercise
  date, the node whose dual cell contains the kink (located by bisection on a cubic interpolant of `C` and the
  exact `E`) receives the cell average of `max(C, E)` (Gauss-Legendre on both sides of the kink); all other nodes
  keep point values (Pooley, Vetzal, Forsyth 2003). This turns the erratic O(dr^2) convergence of point sampling
  into a clean second-order expansion (error ratio 4.00), so Richardson extrapolation is valid. The first, uniform
  scheme without smoothing remains available (`legacyScheme()`); its error changes sign under refinement.
* **Quadrature reference** for the Bermudan option (`CirQuadratureBermudan`, independent of the FD solver):
  backward induction over the exercise dates only. Under the `e_{j+1}`-forward measure
  `C_j(r) = P(e_j, e_{j+1}; r) E^{e_{j+1}}[V_{j+1}(r_{e_{j+1}}) | r]`, and `2 (rho + psi) r_{e_{j+1}}` is noncentral
  chi-square with `4 theta mu / sigma^2` degrees of freedom and noncentrality `2 rho^2 r e^{h tau} / (rho + psi)` -
  the law behind the closed-form zero-bond option, here used for every period. `C_j` is stored as a piecewise
  Chebyshev interpolant on `[0, 1]` (64 panels x 16 nodes); `V_{j+1} = max(C_{j+1}, E_{j+1})` is integrated with
  adaptive Gauss-Kronrod 7/15 between break points (kinks `C = E` and `B = K`, panel boundaries), so each piece is
  smooth. Density: Poisson mixture with Loader's saddle-point terms (validated against 40-digit mpmath values to
  2e-14 relative). Zero bonds in the `LOG1P_RATIO` form. Accuracy ~1e-14 (settings variation), confirmed by an
  independent 30-digit mpmath implementation (`validation/bermudan_quadrature.py`: Bessel-form density, tanh-sinh
  quadrature, Chebyshev-Lobatto panels) and by the extrapolated FD price (RESULTS 6a).

Reference values are generated by `validation/reference_values.py` (output committed in
`validation/reference_values.out`), never by the Java implementation itself.

## 10a. Method comparison against the exact solution (`methodstudy`)

Goal: compare continuation estimators (quantile clusters, least squares) without the confounders of a price
comparison (valuation noise, Euler bias, a reference with its own error, a single training seed).

* **Exact training paths** (`ExactPathGenerator`): the short rate is drawn only at the exercise dates from its exact
  transition law under the forward measure of the last exercise date `S = e_m` (numeraire `P(t, S)`): scaled
  noncentral chi-square with `B(S - T)` in the scale, sampled as Poisson mixture of gamma variables
  (Marsaglia-Tsang). Discounting uses the numeraire, `D(e_j, tau) = P(e_j, S; r_j) / P(tau, S; r_tau)`, stored as
  `cumulativeIntegral = ln(P(e_j, S)/P(0, S))`, so the unchanged induction estimates the same `C_j(r)` as under the
  risk-neutral measure; only the distribution of the training states is the `Q^S` law. The transition law was
  checked with mpmath (mass 1, tower property `P(t,S) E^S[P(T,U)/P(T,S)] = P(t,U)` to 1e-31) and by Monte Carlo
  (zero bonds and Europeans within 4 SE, unit test). For contrast, the same study runs on Euler paths
  (12 steps/year, risk-neutral, rectangle discounting).
* **Exact policy valuation** (`CirQuadratureBermudan.valuePolicy`): a trained rule "exercise iff `E > 0` and
  `E > C^(r)`" is a deterministic function of the rate at each date. Its value in the continuous-time model follows
  from the same quadrature backward induction as the price, with `V_j = E_j` where the rule exercises and the
  rule's own continuation value elsewhere; the rule's switch points (found on a 20 000-point scan plus bisection)
  are integration break points. Settings 32 panels x 12 nodes (price error ~1e-14). The loss
  `L = V* - V(rule) >= 0` is the price a holder following the rule gives away; a threshold at the exact boundary has
  loss < 1e-11 and shifted thresholds lose quadratically in the shift (unit tests).
* **Design:** for each path source, path count `N` and replication (independent training seed
  `seed + 7919 rep`), all estimator configurations are trained on the same paths (paired). Default configurations:
  clusters K = 5, 10, 20, 40 (min. 20 paths, estimation set ALL); LSM degree 2 on ALL paths (the lab's standard),
  degree 2, 3, 4 on in-the-money paths. Paired differences are reported against LSM degree 2 / ALL.
* **Metrics:** policy loss (mean, SE, median, 90% quantile); boundary error `r^_j - r*_j` (switch point closest to
  `r*_j`) and the share of dates where the rule is not a single threshold; continuation error at the exact boundary
  `C^_j(r*_j) - C*_j(r*_j)` (the decision-relevant error) and RMS over the training states (against the optimal
  `C*_j`; the estimators target the continuation under their own later decisions, so this is a diagnostic);
  in-sample price minus `V*`.

## 11. Performance measurement

* Console timings: `System.nanoTime()` single runs without warm-up plus bytes allocated by the calling thread
  (`ThreadMXBean`); orientation only. Logging/CSV are outside the measured blocks; metrics are computed after timing.
* JMH (`lab.cir.bench.PrecisionBenchmarks`): separate areas (analytic zero bond, simulation incl. normal generation,
  exercise values, sorting/partition, backward induction, end-to-end) for precision 16/34/64; warm-up 3 x 1 s,
  5 x 1 s measurements, 1 fork, results returned to JMH (no dead-code elimination), `-prof gc` for allocation per
  operation. Inputs for isolated areas are prepared in `@Setup`; end-to-end includes random number generation.
* Portfolio: independent valuation (own paths per contract, `SE = sqrt(sum SE_c^2)`) vs shared market simulation
  (one path set on the union of exercise dates, per-contract exercise values/policies, portfolio SE from path-wise
  aggregated cashflows). Peak heap via memory pool peaks after `System.gc()` (indicative).

## 12. Known limits

* The Euler + rectangle discretisation has O(dt) bias; quantile clusters and the quadratic LSM basis have
  approximation bias; both are orders of magnitude larger than the arithmetic differences (RESULTS).
* Random inputs are binary64 in both arithmetics; BigDecimal does not remove the randomness' representation error.
* The analytic option references are double based; the FD reference has its own discretisation error (~1e-8 at
  8000 x 3200 with the default scheme, ~1e-11 after extrapolation); the quadrature reference is accurate to ~1e-14
  for the standard case but is so far only validated for `4 theta mu / sigma^2 = 8` (smooth density at r = 0).
* Timings in console/CSV are single-run orientation values; use JMH for statements about performance.
* No Java 8 compatibility module, no Milstein/exact scheme, no trapezoidal discounting control yet (see PLAN).

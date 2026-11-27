# Results (measured, 2026-11-27)

All numbers below come from runs in this repository; experiment ids refer to `results/<id>/` (CSV +
`run_metadata.properties`) and the console logs in `results/logs/`. `results/` is not under version control;
`results/run_batch1.sh` and `results/run_batch2.sh` reproduce the runs (deterministic seeds: prices and
differences are reproducible bit for bit, timings are not).

Environment: OpenJDK 21.0.12 (Ubuntu), WSL2 Linux 6.6 x86_64, 8 logical CPUs, default G1 GC, `-Xmx4g`;
big-math 2.3.2, commons-math3 3.6.1. **Timings in this section are orientation values** (single runs,
`System.nanoTime`) unless marked as JMH.

Standard case (spec 19.1): `theta = 0.5, mu = 0.04, sigma = 0.10, r0 = 0.03` (Feller: 0.04 >= 0.01),
5y annual 4% bond on 100, strike 100, exercise dates 1, 2, 3, 4; 12 Euler steps/year; training seed 42,
valuation seed 4242; 10 quantile clusters (min. 20 paths); LSM basis `1, r, r^2`; estimation set ALL;
BigDecimal 34 digits HALF_EVEN; zero-bond form `LOG_DIFFERENCE`.

## 1. Headline

| quantity (standard case unless stated) | size |
|---|---|
| precision effect `V_dec - V_dbl`, out-of-sample, cluster (1 000 paths) | **6.8e-14** (3.5e-14 relative, 1.5e-12 SE) |
| same with 20 000 + 20 000 paths, 40 clusters (research run) | 7.9e-14 (8.0e-12 SE) |
| Monte Carlo standard error (1 000 / 20 000 paths) | 4.4e-2 / 9.8e-3 |
| training-seed variation of the out-of-sample price (8 seeds, 1 000 paths) | 1.908 ... 1.917 (spread ~9e-3) |
| method effect `V_LSM - V_Cluster`, out-of-sample (20 000 paths, one seed, Euler model) | -1.7e-3 (2.3% of paths exercise on a different date); **not confirmed** by the exact evaluation, section 6b |
| policy loss `V* - V(rule)`, exact valuation, 16 000 exact training paths, mean over 40 seeds | LSM degree 2 / ITM 3.6e-5, best clusters (K = 10) 2.5e-4, section 6b |
| Euler/rectangle bias at 12 steps/year vs deterministic reference (20 000 paths, out-of-sample) | +2.3e-2 (4 steps: +7.6e-2; 252 steps: +2.1e-3) |
| exercise decisions that differ between double and BigDecimal (34 digits), any run up to 20 000 paths | **0** |
| cluster switches in Mode B (34 digits), any run | **0** |
| runtime factor BigDecimal / double, end to end (JMH, 34 digits) | ~230 (16 digits: ~130, 64 digits: ~340); allocation ~1 300x |

`|Delta_method| / |Delta_precision| = 2.1e10` in the research run (`20261127-093635-889-compare-research`).
The ordering of error sources found here is

```
precision (1e-13)  <<  MC noise (1e-2)  ~  Euler/discretisation bias (1e-2 at 12 steps/year)
                   ~  policy training noise (1e-2)  >  method choice cluster vs LSM (policy loss 1e-4..1e-3, section 6b)
```

for this contract; i.e. hypotheses 1, 4, 5 (relative to precision) and 8 are supported, hypotheses 2 and 3 are not
observable at 34 digits (but demonstrably real at <= 6 digits, section 5).

## 2. Where the double-vs-BigDecimal difference comes from (`20261127-092857-608-compare-standard`)

State differences at the exercise dates (training paths, BigDecimal - double):

| date | max abs dr | max ULP of r | max abs dD (discount) | max rel dD | max abs dB (bond value) | max abs dE |
|---|---|---|---|---|---|---|
| 1 | 3.2e-17 | 6.9 | 6.5e-17 | 6.7e-17 | 1.7e-13 | 1.7e-13 |
| 2 | 3.8e-17 | 6.4 | 9.8e-17 | 1.1e-16 | 7.0e-14 | 7.0e-14 |
| 3 | 5.2e-17 | 6.1 | 1.3e-16 | 1.4e-16 | 2.4e-14 | 2.4e-14 |
| 4 | 4.5e-17 | 7.8 | 1.5e-16 | 2.1e-16 | 6.9e-14 | 6.9e-14 |

* **Simulation:** the Euler short rate drifts only ~1 ULP on average (max ~8 ULP after 48 steps); no growth beyond
  that along the time grid (mean ULP 0.3 at t=0 [input rounding of r0], ~1.0-1.3 for t >= 1).
* **Discounting:** relative error of `exp(-I)` ~1e-16.
* **Exercise values dominate:** `max |dE| = 1.7e-13` on values ~1-5, i.e. ~1e-15 relative *of the bond (~100)*.
  The cause is cancellation in the textbook form of `ln A(tau)` of the analytic zero bond (relative error of P up to
  ~2e-15 for tau = 4, measured against mpmath), multiplied by the bond notional and then exposed by `B - K`.
  With the algebraically equivalent `LOG1P_RATIO` form, `Delta V` falls from 6.8e-14 to **1.4e-14**
  (`20261127-093335-978-sensitivity-zerobond.form-standard`) - a pure double-side fix of 80% of the gap.
* **Cluster continuation values:** max `|C_dec - C_dbl| = 5.2e-14`, i.e. 7e-13 of the cluster's own standard error
  (`s/sqrt(n)`): numerically irrelevant relative to the statistics inside a cluster.
* **Regression (LSM, `1, r, r^2`):** coefficients differ by at most 3.9e-14 relative.
* **Frozen double policy with BigDecimal cashflows** gives the same paired difference as the full comparison
  (6.77e-14): with zero decision differences, the whole price gap is cashflow arithmetic.
* **Double aggregation error:** `Delta V` (6.77e-14 out-of-sample) differs from the exact paired mean (6.77e-14 at
  1 000 paths, but 7.07e-14 vs 7.85e-14 at 20 000 paths) - the naive double summation of the path values contributes
  ~8e-15 at 20 000 paths and grows with the path count (paths series). Kahan summation of the *cluster means* made no
  difference at all (clusters have only ~100-500 terms); the final mean is the place where compensated or pairwise
  summation would act.

Paired statistics (12.11), standard case, cluster, out-of-sample: mean difference 6.77e-14,
`SE_d = 2.5e-15`, 95% CI [6.28e-14, 7.25e-14], max |d_p| = 1.8e-13, 827 of 1 000 paths with d != 0 (all others are
paths without exercise), 0 paths with a different exercise date. The difference is statistically "significant"
(CI excludes 0) - it is systematic rounding, not noise - but 12 orders of magnitude below the price SE.

## 3. Comparison modes (A vs B) and decisions

In every run with >= 7 digits, Mode A and Mode B produced identical policies: no state crossed a quantile
boundary (Mode B cluster switches = 0), no exercise decision differed (in-sample or out-of-sample), no target
changed through propagation (12.8: 0 paths with a different later exercise index at every date).
This holds for paths 500-20 000, cluster counts 5-50, steps 4-84/year, strikes 92-108, volatilities 0.05-0.40
(incl. Feller violations), ITM/ALL selection, both LSM bases, training seeds 1-8 and 1-8 threads.

Why: for a flip, a path's state must lie within ~|dr| ~ 5e-17 of a cluster boundary, or `|E - C|` must be below
~|dE| ~ 1e-13. With continuous state distributions the expected number of such paths is
~(number of paths) x (density) x (difference): e.g. 2e4 paths x O(1) per unit of `E - C` x 1e-13 ~ 2e-9.

Approximate exercise boundaries (12.9) in the short rate, out-of-sample, 20 000 paths:
`e1: [0.029702, 0.029704]` (cluster) vs `[0.029031, 0.029033]` (LSM); `e2: 0.03116` vs `0.03101`;
`e3: 0.03424` vs `0.03404`; `e4: 0.03907` (no continuation). Double and BigDecimal boundaries are identical to all
printed digits; the two *methods* differ by 2-7e-4 in the rate. The brackets are narrow because exercise regions
turned out to be intervals `r < r*` for both methods (overlap count 0).

## 4. Precision series (`20261127-093250-339-sensitivity-decimal.precision-standard`, 1 000 paths)

BigDecimal price (Mode B, out-of-sample) minus the 100-digit run of the same algorithm:

| digits | 8 | 10 | 12 | 16 | 20 | 34 | 50 | 64 |
|---|---|---|---|---|---|---|---|---|
| cluster | -1.7e-5 | 4.3e-8 | 1.1e-9 | -1.8e-13 | 2.8e-17 | 1.5e-31 | 2.9e-48 | -5.4e-63 |

* 34 digits reproduce the 100-digit reference to ~31 digits; 64 vs 100 agree to ~63 digits (stable reference).
* double (`Delta V` vs the 34-digit run: 6.8e-14) sits between the 16- and 20-digit decimal runs, as expected for
  binary64 (~15.95 digits) with a few ULP of error accumulation.
* Rounding modes at 16 digits (`...decimal.roundingMode...`): HALF_EVEN -1.2e-13, HALF_UP -4.9e-14,
  DOWN -1.2e-13, CEILING +4.7e-13; directed rounding (CEILING) produces a visibly biased drift.

## 5. Constructed stress test: discrete amplification at very low precision (`...stress-lowprec`, 4 000 paths)

To make hypotheses 2/3 observable, the decimal precision was lowered artificially (these runs are *not*
realistic configurations):

| digits | `Delta V` (cluster) | `|Delta V|/SE` | decision differences (in + out) | cluster switches (Mode B) | LSM decision differences |
|---|---|---|---|---|---|
| 7 | -1.7e-4 | 0.008 | 0 | 0 | 0 |
| 6 | 2.1e-3 | 0.09 | 3 | 0 | 10 |
| 5 | -9.4e-3 | 0.42 | 155 | 4 | 168 |
| 4 | 9.8e-2 | 4.4 | 1 094 | 32 | 1 109 |

Discrete effects (flipped decisions, cluster switches, propagation) appear only once the arithmetic error of the
exercise values reaches ~1e-4..1e-3 absolute (relative ~1e-6). They exist, they propagate, and at 4 digits they
dominate the MC error - but double precision is ~9 orders of magnitude away from that regime.

## 6. Error budget against independent references

**European control cases** (`20261127-093516-998-european-european50k`, 50 000 valuation paths, seed 777):

| case | reference | FD | MC double | SE | z | dec - dbl |
|---|---|---|---|---|---|---|
| ZBC, unit zero bond S=5, T=1, X=0.85 | 0.015910 (closed form) | 0.015910 | 0.015986 | 6.5e-5 | 1.2 | 1.0e-15 |
| coupon bond call T=1, K=100 | 1.423263 (Jamshidian) | 1.423272 | 1.433907 | 6.6e-3 | 1.6 | 1.0e-13 |
| T=2 | 1.204283 | 1.204290 | 1.208885 | 6.1e-3 | 0.8 | 3.5e-14 |
| T=3 | 0.921652 | 0.921649 | 0.934899 | 4.9e-3 | 2.7 | -8.5e-15 |
| T=4 | 0.546975 | 0.546978 | 0.553438 | 3.0e-3 | 2.2 | -2.9e-14 |

Euler step series for the European at T=4 (coupled increments, same 252/year Brownian path):
MC - Jamshidian = 1.6e-2 (4/year), 8.0e-3 (12), 5.1e-3 (36), 4.4e-3 (84), 4.0e-3 (252, z = 1.4).
The Euler + rectangle-rule bias at 12 steps/year is thus ~4e-3..8e-3 for T = 4 and positive; the z-values around
2-3 at 12 steps are consistent with this bias, not with a pricing error (the FD column was produced with the first,
uniform FD scheme; FD and Jamshidian agree to ~1e-5 there, to ~1e-7 with the current default scheme, section 6a).

**Bermudan vs deterministic reference** (`20261127-092958-036-reference-standard`): the Monte Carlo prices below
were compared with the first FD reference 1.863362 (uniform grid, 4000 x 1600). The accurate value is
1.8633604211 (section 6a); the differences in the table change by 2e-6, below the digits shown.

| steps/year | cluster out-of-sample - reference | LSM out-of-sample - reference | SE |
|---|---|---|---|
| 4 | +7.6e-2 | +7.8e-2 | 1.0e-2 |
| 12 | +2.3e-2 | +2.4e-2 | 9.8e-3 |
| 36 | +9.3e-3 | +8.7e-3 | 9.7e-3 |
| 84 | +4.2e-3 | +4.1e-3 | 9.6e-3 |
| 252 | +2.1e-3 | +2.1e-3 | 9.6e-3 |

(20 000 training + 20 000 valuation paths, coupled increments.) For the same setup at 12 steps/year,
`V_dec - V_dbl` = 7.8e-14 (cluster) / 7.4e-14 (LSM). The out-of-sample price is a lower bound only for the
*discretised* model; the positive Euler bias of the exercise values dominates the (negative) policy sub-optimality.
The European price at e1 (1.4233) and all other Europeans are clearly below the Bermudan (1.8634).

## 6a. Accuracy of the deterministic Bermudan references (`20261127-150610-637-fdstudy-standard`)

Two independent deterministic methods for the continuous-time Bermudan price (log: `results/logs/fdstudy.txt`):

* **Quadrature** (`CirQuadratureBermudan`): backward induction over the four exercise dates only; continuation value
  `C_j(r) = P(e_j, e_{j+1}; r) E^{e_{j+1}}[V_{j+1}]` under the forward measure, where the transition law is exactly a
  scaled noncentral chi-square; piecewise Chebyshev interpolant for `C_j`, adaptive Gauss-Kronrod between the kinks.
  No time-stepping, no Euler bias.
* **Finite differences** (`CirFiniteDifferenceBermudan`), Crank-Nicolson + Rannacher, in four variants.

| quadrature settings (rMax, panels x nodes, tol) | price |
|---|---|
| 1.0, 32 x 12, 1e-11 (1.2 s) | 1.863360421058373 |
| 1.0, 64 x 16, 1e-13 (standard, 5.8 s) | 1.863360421058385 |
| 1.0, 128 x 20, 1e-13 (19 s) | 1.863360421058386 |
| 1.5, 96 x 16, 1e-13 | 1.863360421058385 |
| mpmath, 30 digits, independent implementation (`validation/bermudan_quadrature.py`, 24 x 15 Lobatto, 47 min) | 1.863360421058401499 |

The double quadrature agrees with the 30-digit mpmath value to **1.6e-14** (Bermudan) and to 2e-15 ... 3e-14 for the
four Europeans; the exercise boundaries agree to < 4e-16. Against the double Jamshidian closed form the European
differences are 3e-14 ... 9e-14 - the closed form (double noncentral chi-square CDF) is the less accurate of the two.
The mpmath value is limited by its own interpolation (24 panels x 15 nodes), which was not refined further
(runtime); its agreement with the double code to 1.6e-14 bounds both.

FD under joint refinement (rate intervals x steps/year), error against the quadrature price:

| grid | uniform (first scheme) | uniform + kink smoothing | sinh grid | **sinh + smoothing (default)** |
|---|---|---|---|---|
| 250 x 100 | +6.6e-4 | +1.1e-3 | -2.3e-5 | +9.3e-6 |
| 500 x 200 | +1.3e-4 | +2.6e-4 | -4.7e-6 | +2.3e-6 |
| 1000 x 400 | +2.3e-5 | +6.6e-5 | -2.8e-6 | +5.7e-7 |
| 2000 x 800 | **-3.6e-6** | +1.7e-5 | -3.3e-7 | +1.4e-7 |
| 4000 x 1600 | +1.9e-6 | +4.1e-6 | -5.1e-8 | +3.6e-8 |
| 8000 x 3200 | +4.0e-7 | +1.0e-6 | -4.6e-8 | +8.9e-9 |
| error ratio per doubling | 4.9, 5.9, -6.4, -1.9, 4.9 | 4.05, 4.00, 4.01, 4.00, 4.00 | 4.8, 1.7, 8.4, 6.6, 1.1 | 4.04, 4.07, 4.01, 3.97, 4.00 |
| Richardson (finest pair) - quadrature | -1.1e-7 | -3.5e-10 | -4.4e-8 | **+6.8e-12** |

* **The first FD scheme did not converge monotonically**: its error changes sign between 1000 and 4000 intervals,
  so the change between grids (5.5e-6 at the finest pair, reported earlier) was not an error bound; its actual
  error at 4000 x 1600 was +1.9e-6. Cause: the kink of `max(C, E)` (and of the final payoff) lies at an arbitrary
  position inside a grid cell; point sampling gives an O(dr^2) error whose constant depends on that position.
* **Cell averaging at the kink** (the node whose dual cell contains the kink receives the cell average of
  `max(C, E)`) restores a clean ratio of 4.00, i.e. a smooth error expansion; Richardson extrapolation then becomes
  valid. On the uniform grid this costs accuracy at fixed size (the error is larger but regular).
* The **sinh grid** concentrates nodes near `r0` and the exercise boundaries (0.029 ... 0.039): 30-70x smaller
  errors, but without smoothing the convergence is again erratic.
* With both, FD (8000 x 3200, Richardson) and quadrature agree to **7e-12** (to the mpmath value likewise 7e-12): two
  independent methods confirm the Bermudan price to about 11 significant digits.

Separate refinement (default scheme): rate grid only (6400 steps/year) gives errors +2.1e-5 ... +1.8e-8 with
ratios 4.0; time grid only (16000 intervals) gives -1.9e-4 ... -4.1e-8 with ratios 4.0 (the last ratios 4.5/4.4
contain the residual error of the fixed dimension). **Space and time errors have opposite signs**, which is why the
joint refinement error is smaller than either. The `rMax` series (0.2 ... 1.0 at 4000 x 1600) changes the price by
5e-8, but this reflects the coarser resolution near `r0` for a longer sinh grid, not truncation: the extrapolated
value at `rMax = 0.5` already agrees with the quadrature to 7e-12.

Exercise boundaries `r*_j` (exercise below): 0.029278343949, 0.031428653031, 0.034423795917 at e1..e3 and
0.039067630129 (bond value = strike) at e4; FD (8000 x 3200) agrees to 6e-11.

**Consequence for method comparisons.** The deterministic reference error (~1e-11) is now 8 orders of magnitude
below the Monte Carlo effects of section 6, and the quadrature also delivers the exact continuation functions
`C_j(r)` and boundaries, so a policy estimated by clusters or LSM can be valued without simulation noise
(next step, see PLAN).

## 6b. Cluster vs LSM against the exact solution (`20261127-163245-140-methodstudy-standard`)

Design (METHODOLOGY 10a): 2 training-path sources (exact simulation at the exercise dates; Euler 12 steps/year)
x 7 path counts x 40 independent training seeds x 8 estimator configurations, all configurations trained on the same
paths; every trained rule valued exactly by quadrature (4 480 rules, 12 min on 8 threads; log
`results/logs/methodstudy.txt`, raw and summary CSV in the run directory). `V* = 1.863360421058`, exact boundaries
`r* = 0.029278, 0.031429, 0.034424` (e1..e3).

**Policy loss `V* - V(rule)`, exact training paths**, mean over 40 seeds (standard error in parentheses):

| N | CL5 | CL10 | CL20 | CL40 | LSM2 ALL | LSM2 ITM | LSM3 ITM | LSM4 ITM |
|---|---|---|---|---|---|---|---|---|
| 250 | 8.1e-3 (1.1e-3) | 1.0e-2 (1.4e-3) | 1.4e-2 (1.6e-3) | = CL20 * | **1.5e-3** (2.5e-4) | 2.4e-3 (4.2e-4) | 3.6e-3 (5.2e-4) | 4.9e-3 (6.6e-4) |
| 1 000 | 3.2e-3 (4.2e-4) | 2.6e-3 (3.5e-4) | 4.7e-3 (5.5e-4) | 8.0e-3 (6.9e-4) | **4.7e-4** (7.4e-5) | 5.2e-4 (8.3e-5) | 7.2e-4 (9.6e-5) | 8.2e-4 (1.1e-4) |
| 4 000 | 2.4e-3 (2.7e-4) | 9.2e-4 (2.0e-4) | 1.6e-3 (2.2e-4) | 2.9e-3 (3.1e-4) | 2.4e-4 (3.1e-5) | **1.4e-4** (2.3e-5) | 2.1e-4 (3.8e-5) | 2.3e-4 (3.9e-5) |
| 16 000 | 1.9e-3 (1.0e-4) | 2.5e-4 (2.9e-5) | 3.9e-4 (4.8e-5) | 7.7e-4 (9.4e-5) | 1.5e-4 (9.7e-6) | **3.6e-5** (6.7e-6) | 4.2e-5 (5.9e-6) | 5.0e-5 (5.5e-6) |

\* with at least 20 paths per cluster, 250 paths allow at most 12 clusters. Full table (also N = 500, 2 000, 8 000,
medians, 90% quantiles, paired differences): log file.

**Diagnostics at N = 16 000 (exact paths):**

| | CL5 | CL10 | CL20 | CL40 | LSM2 ALL | LSM2 ITM | LSM3 ITM | LSM4 ITM |
|---|---|---|---|---|---|---|---|---|
| mean abs boundary error `|r^ - r*|` | 3.6e-4 | 1.8e-4 | 2.1e-4 | 2.1e-4 | 2.3e-4 | **8.3e-5** | 9.8e-5 | 1.1e-4 |
| dates with more than one switch point | 33% | 28% | 38% | 47% | 0% | 0% | 0% | 0% |
| RMS continuation error at `r*` | 6.9e-2 | 3.8e-2 | 4.5e-2 | 5.2e-2 | 2.7e-2 | **1.2e-2** | 1.4e-2 | 1.5e-2 |
| paired loss difference vs LSM2 ALL | +1.8e-3 | +9.8e-5 | +2.4e-4 | +6.2e-4 | - | **-1.1e-4** (95% of seeds better) | -1.1e-4 | -9.9e-5 |

**Euler training paths** (12 steps/year), N = 16 000: CL5 3.0e-3, CL10 6.4e-4, CL20 3.7e-4, CL40 6.4e-4,
LSM2 ALL 3.3e-4, LSM2 ITM **7.1e-5**, LSM3 ITM 1.0e-4, LSM4 ITM 1.1e-4; in-sample price minus `V*` +2.6e-2 for all
methods (the Euler bias of section 6, independent of the estimator).

Findings:

1. **The best LSM configuration beats the best cluster configuration at every path count**, in the mean, the median
   and the 90% quantile, by a factor of 5-10 in the mean loss (16 000 exact paths: 2.5e-4 vs 3.6e-5).
   The cause is visible in the diagnostics: a cluster estimate is constant over the cluster, so near the boundary it
   averages a steep continuation function (RMS continuation error at `r*` 3x larger, best vs best). That shifts the
   switch point, and noisy cluster means produce non-monotone exercise sets (28-47% of dates have several switch
   points; LSM rules are thresholds).
2. **Fixed K does not converge.** CL5 stays at ~2e-3 from 4 000 paths on (pure approximation bias); the best K grows
   with N (K = 5 at 250-500 paths, K = 10 from 1 000; K = 20/40 are worse at every N tested). A cluster count rule
   must therefore depend on N - and even the best K is behind LSM.
3. **LSM estimation set:** ALL is better for small samples (N <= 500: more data, less variance; significant in the
   paired differences), the two are equal at 1 000-2 000 paths, ITM is better from 4 000 paths. With ALL the loss saturates at ~1.5e-4 (the quadratic is fitted over the whole state range,
   which biases it at the boundary); with ITM it keeps falling to 3.6e-5. Degrees 3 and 4 do not help here: the
   continuation function in the in-the-money range is close to quadratic, and extra terms only add variance.
4. **Euler training paths** roughly double the loss of the best rules (7.1e-5 vs 3.6e-5 at 16 000) - the rule is
   trained on a slightly wrong model - but do not change the ranking.
5. **Loss vs boundary error:** losses are roughly quadratic in the boundary error (unit test for shifted thresholds;
   here e.g. LSM2 ITM 8e-5 boundary error / 3.6e-5 loss, CL10 1.8e-4 / 2.5e-4). This is why price comparisons are
   insensitive: a boundary that is visibly wrong costs little price.
6. **Earlier single-seed price comparison revisited.** The out-of-sample comparison of section 1 (one seed, 20 000
   paths, K = 40 vs LSM2 ALL, valued by Monte Carlo in the 12-step Euler model) found the cluster policy *higher*
   by 1.7e-3. The exact evaluation over 40 seeds shows the opposite on average (Euler paths, 16 000: CL40 loses
   3.1e-4 more than LSM2 ALL, cluster better in only 23% of seeds). The single-seed result was within the seed
   spread, and it measured the value in the discretised model rather than in the continuous model.
7. All effects here (1e-5 ... 1e-2) are 8-11 orders of magnitude above the double-vs-BigDecimal effect (1e-13).

## 7. Further sensitivity results (out-of-sample, cluster, 1 000 paths unless stated)

| factor | values | precision effect `Delta V` | remark |
|---|---|---|---|
| paths (train = valuation) | 500 ... 16 000 | 7.0e-14 ... 7.9e-14 | `|Delta V|/SE` grows like sqrt(n) but stays < 1e-11 |
| clusters (4 000 paths) | 5, 10, 20, 50 | ~7e-14 | prices 1.9263, 1.9259, 1.9260, 1.9267 (method effect of K ~ 1e-3 < SE) |
| strike | 92, 96, 100, 104, 108 | 1.0e-13, 8.3e-14, 6.8e-14, 1.3e-14, 0 | K=108 > bond value at r=0 for all paths: price 0 |
| volatility | 0.05 ... 0.40 | -1.3e-13 ... 2.4e-15 | sigma 0.30/0.40 violate Feller (ALLOW_VIOLATION) and produce negative Euler states (identical counts in both arithmetics is a unit test); even at sigma = 0.10 (Feller satisfied) the research run had 3 negative internal states on 20 000 paths, identically in both arithmetics |
| estimation set | ALL / ITM | 6.8e-14 / 6.7e-14 | prices 1.9166 / 1.9203 (cluster), 1.9162 / 1.9145 (LSM) |
| LSM basis | r / B | 6.7e-14 / 6.7e-14 | the badly conditioned `1, B, B^2` basis did not increase the precision effect |
| Euler steps (coupled) | 4, 12, 36, 84 | 6.1e-14 ... 6.4e-14 | precision effect independent of the grid |
| threads (8 000 paths) | 1, 2, 4, 8 | identical results | BigDecimal 11.3 s -> 5.9 s (only simulation is parallel) |

## 8. Cost

Orientation timings, standard case, 1 000 paths, single thread (`20261127-092857-608-compare-standard`):

| step | double | BigDecimal | factor | allocated (dbl / dec) |
|---|---|---|---|---|
| paths: simulation + exercise values (1 000 x 48 steps) | 8 ms | 1 070 ms | ~130 | 0.6 MB / 884 MB |
| cluster backward induction | 4 ms | 88-164 ms | 20-37 | 0.2 MB / 120 MB |
| LSM backward induction | 3 ms | 113-123 ms | 33-36 | 0.2 MB / 157 MB |
| out-of-sample valuation | 0.6-1.6 ms | 26-31 ms | 16-48 | 0.1 MB / 35-40 MB |

Cost grows with the precision (precision series, total BigDecimal time for 1 000 + 1 000 paths): 1.2 s (16 digits),
1.4 s (34), 2.4 s (64), 3.4 s (100).

**Portfolio** (`...portfolio-10...`, `...portfolio-100...`, 1 000 paths, in-sample, clusters):

| contracts | strategy | double | BigDecimal | contracts/s (dbl / dec) | allocated (dbl / dec) | peak heap (dbl / dec) |
|---|---|---|---|---|---|---|
| 100 | independent | 1.5 s | 49.9 s | 67 / 2.0 | 0.3 GB / 72 GB | 74 MB / 275 MB |
| 100 | shared simulation | 0.07 s | 22.8 s | 1 530 / 4.4 | 20 MB / 37 GB | 20 MB / 284 MB |

Shared simulation pays off strongly in double (23x) but only 2.2x in BigDecimal, because per-contract work
(BigDecimal bond values with one `exp` per cashflow and path, induction) dominates. Portfolio value differences:
dec - dbl = 9.9e-12 (independent) / 9.2e-12 (shared) on a portfolio value of ~589. The shared-simulation portfolio SE
(5.97) is much larger than the independent one (1.05) because all contracts share paths and are strongly positively
correlated - the path-wise aggregated SE accounts for that; adding individual SEs would not.

JMH figures: section 9.

## 9. JMH (isolated measurement areas)

`java -cp target/cir-numerical-precision-lab.jar org.openjdk.jmh.Main lab.cir.bench -prof gc`, JMH 1.37,
average time, 3 x 1 s warm-up, 5 x 1 s measurement, 1 fork, `-Xmx2g`, nothing else running; 1 000 paths,
48 steps, 4 exercise dates; raw data in `results/jmh-results.csv`, log in `results/logs/jmh.txt`.

| area | double ms/op | BigDecimal 16 | 34 | 64 digits | factor at 34 | allocated B/op double | BigDecimal 34 |
|---|---|---|---|---|---|---|---|
| analytic zero bond (coefficients + price, one call) | 7.1e-5 | 0.139 | 0.254 | 0.401 | ~3 600 | ~0 (escape analysis) | 4.5e5 |
| simulation incl. normal generation | 2.14 | 173 | 327 | 439 | 153 | 5.8e5 | 4.9e8 |
| exercise values (bond value on all paths/dates) | 0.141 | 157 | 262 | 470 | ~1 900 | ~1 | 4.4e8 |
| sorting + quantile partition | 0.42 | 0.68 | 0.91 | 0.84 | 2 | 1.2e5 | 1.2e5 |
| backward induction (clusters) | 0.081 | 51 | 80 | 171 | ~1 000 | 1.2e5 | 1.3e8 |
| end to end (normals ... price) | 3.16 | 414 | 729 | 1 086 | 231 | 8.2e5 | 1.05e9 |

Score errors (99.9% CI) are ~2-15% for most areas but large for the allocation-heavy BigDecimal induction and
clustering (up to ~80% at 16 digits) and end-to-end at 64 digits (~36%): GC noise with only 5 iterations and 1 fork.
The factors are therefore order-of-magnitude statements; more forks/iterations are needed for firm numbers.
Observations:

* End to end, BigDecimal costs a factor ~130 (16 digits), ~230 (34) and ~340 (64) and allocates ~1 300x more.
* The most expensive single area is the exercise value: one big-math `exp` per cashflow and path
  (~25 us per `exp` at 34 digits vs ~15 ns in double). The simulation (one `sqrt` per step) is second.
* Sorting BigDecimal states is only 2x slower than doubles (dominated by boxing in the comparator sort in both).
* Caveat: the zero-bond micro benchmark uses constant arguments; for double, the JIT removes the allocation of the
  coefficient record (escape analysis) - the ~70 ns are plausible for `expm1`, 2 `log`, `exp`, but the micro number
  should not be over-interpreted.

## 10. Conclusions for the decision "when is BigDecimal worth it?"

1. For valuation of this product class with double-precision random inputs, **BigDecimal does not change any decision
   and changes prices by ~1e-13**, 11-12 orders of magnitude below MC noise, Euler bias and policy/method effects,
   at a cost factor of ~230 in time and ~1 300x in allocation (JMH, 34 digits).
2. The largest double error is not in the simulation but in a **numerically unfavourable closed-form expression**
   (cancellation in `ln A`). Fixing the formula in double removes most of the gap at zero cost
   (hypothesis 7 confirmed for this case). Compensated summation of the final mean is the next candidate.
3. Discrete amplification (flipped decisions, cluster switches, propagation) is real but needs relative errors of
   ~1e-6 or worse; in double it is not observable even with 20 000 paths and 40 clusters.
4. Better spent effort: finer Euler grids or a better scheme (bias 2e-2 at 12 steps/year), more paths, out-of-sample
   policies and method checks against a deterministic reference (quadrature / FD, section 6a), the choice of the
   continuation estimator (the best LSM rule loses 5-10x less than the best quantile-cluster rule, section 6b),
   shared simulations for portfolios.
5. BigDecimal remains valuable as an **algorithm reference** (34 digits reproduce 100 digits to 31 digits) for
   locating double error sources - as done here - not as a production arithmetic.

These conclusions are specific to the tested parameter ranges (section 7) and to binary64 random inputs; they are
empirical, not proofs.

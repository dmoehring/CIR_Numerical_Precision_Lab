# Development plan and status

The specification asked for a plan first and phase-wise approval; on request the phases 0-5 plus first versions
of phase 6 and parts of phase 7 were implemented in one go. This document records the plan, what is done, the
acceptance criteria, the conventions that are fixed, the decisions still open, and the risks.

## Phases

| phase | content | status |
|---|---|---|
| 0 | requirements, mathematics, architecture, validation strategy | done (`docs/METHODOLOGY.md`) |
| 1 | CIR model (Feller modes), cashflows, double/BigDecimal zero bonds, common random numbers | done, tested |
| 2 | full-truncation Euler, rectangle discounting, European single-date control case, analytic references (closed form, Jamshidian), Euler step series | done, tested |
| 3 | Bermudan backward induction with quantile clusters, Mode A/B, detailed metrics, console report, first timings | done, tested |
| 4 | out-of-sample policies, paired statistics, frozen-policy control, sensitivity and precision series, CSV, run metadata | done |
| 5 | Longstaff-Schwartz with Householder QR in both arithmetics (`1, r, r^2` or `1, B, B^2`) | done, tested |
| 6 | JMH benchmarks, portfolio (independent vs shared simulation), allocation/peak-heap figures, path-parallel simulation | first version done; thread-scaling series and larger portfolios open |
| 7 | Jamshidian (done), FD Bermudan reference (done), alternative zero-bond evaluation form (done), Kahan control (done); open: high-precision noncentral chi-square, Milstein/exact scheme, trapezoidal rule, hybrid double-simulation/BigDecimal-aggregation, Java 8/21/25 + StrictMath comparison, HTML report | partly done |

## Acceptance criteria for the MVP (all met)

1. `mvn package` builds with Java 21, all tests deterministic and green (47 tests).
2. Feller violation rejected in `STRICT`, allowed in `ALLOW_VIOLATION`; other restrictions in both modes.
3. Zero bonds in both arithmetics agree with independent 60-digit references within a documented error model.
4. Both arithmetics consume the identical binary64 normals (test reproduces the first Euler step by hand).
5. European single-date MC agrees with Jamshidian within 4 SE (20 000 paths, fixed seed); Bermudan out-of-sample is not
   below any European at an exercise date beyond 4 SE; FD Bermudan >= FD European.
6. Constructed backward-induction cases (tie -> continue, no-exercise -> 0, replacement of later cashflows,
   ITM selection, cluster means) give the hand-computed results in both arithmetics.
7. Console report shows price/SE/time per variant, `Delta V` absolute/relative/over SE, paired CI, method name,
   cluster mode, precision; CSV export with full decimal strings and run metadata.
8. Mode A uses identical cluster path indices (tested); Mode B with a constructed swap shows exactly the expected
   switches near the boundary (tested).
9. Results do not depend on the thread count (tested).

## Fixed conventions (see METHODOLOGY for the reasoning)

* Cashflows at an exercise date are not part of the remaining bond; strike is a money amount.
* Ties `E = C` continue; exercise only with `E > 0`; no later exercise -> target 0.
* Full-truncation Euler, effective rate `x^+` everywhere (no projection of the stored state); left rectangle rule.
* Estimation set `ALL` (configurable `IN_THE_MONEY`); selection independent of earlier decisions.
* Quantile clusters on the effective short rate, sort key `(state, pathIndex)`, sizes `floor` split,
  `minPathsPerCluster` reduces K; new states assigned by right-closed intervals with open ends.
* Mode A: structure (eligible set, clusters, bounds) from the double states; Mode B: own structure.
* Common random numbers: SplitMix64 per path + AS 241; `new BigDecimal(z)`.
* DECIMAL128 default, guard digits 10, big-math for exp/log/sqrt; differences computed in BigDecimal.

## Open decisions (for joint discussion)

1. **Default zero-bond form.** `LOG_DIFFERENCE` (textbook) is the default so that the baseline measures what a
   straightforward implementation does. `LOG1P_RATIO` reduces the double error of the exercise values ~3x. Should
   the stable form become the default (and the textbook form a documented control)?
2. **Default estimation set** for LSM: `ALL` (uniform with clusters) vs `IN_THE_MONEY` (literature standard).
3. **Cluster count rule**, e.g. `K = clamp(sqrt(n)/4, 5, 50)` vs fixed K; the cluster-count series shows the
   bias/variance trade-off.
4. Whether low-precision BigDecimal runs (8-12 digits) should stay in the standard precision series as a
   constructed stress test for discrete amplification (they are clearly marked as artificial).
5. Scaling of the LSM basis (`B` centred/scaled) - currently deliberately unscaled.
6. Portfolio and thread grid for the larger performance study (phase 6) and whether a Java-8-compatible core module
   is wanted for the Java 8/21/25 comparison.

## Where the double-vs-BigDecimal comparison is methodologically weak

* **Random inputs are binary64** in both variants. BigDecimal reproduces the *same* discretised experiment with less
  rounding, it does not create a more accurate Brownian path.
* **The algorithm defines the reference.** A 100-digit run is the exact value of *this* Euler/rectangle/cluster
  algorithm on *these* normals, not the price of the option. Its distance to the true price is dominated by MC,
  discretisation and continuation-estimation errors that are 10-12 orders of magnitude larger.
* **Discrete decisions are not differentiable.** Whether a tiny difference flips a decision is a matter of the
  distance to the next cluster boundary or exercise boundary; with continuous state distributions the probability is
  proportional to the size of the difference (~1e-16 relative) times the density at the boundary. Not observing any
  flip in finite samples is therefore expected and not evidence of robustness in constructed degenerate cases.
* **The double baseline is implementation specific.** The size of `Delta V` depends on formula forms
  (`LOG_DIFFERENCE` vs `LOG1P_RATIO`), summation order and `Math` vs `StrictMath`; "double" is not one fixed number.
* **Timing comparison is not like for like at the micro level:** BigDecimal allocates on every operation; costs of
  big-math functions grow super-linearly with digits; JIT and GC behaviour differ. Only JMH figures should be used
  for conclusions.

## Risks

| risk | mitigation |
|---|---|
| memory: BigDecimal paths x dates x 4 arrays + per-date diagnostic arrays | only exercise-date states are stored; normals generated per path; diagnostics per date; research runs started explicitly |
| runtime of BigDecimal (factor ~100 in simulation) | path-parallel simulation with deterministic results; smaller precision series sizes |
| continuation estimation bias / variance | cluster count series, LSM comparison, out-of-sample valuation, FD reference |
| comparability of arithmetic | line-by-line mirrored operation order, shared normals, Mode A, frozen policy control |
| instability of discrete decisions | Mode B diagnostics (switches, rank distances, |E - C| of flips), low-precision stress series |
| validity of references | independent Python/mpmath/SciPy values, documented validated ranges, FD convergence |

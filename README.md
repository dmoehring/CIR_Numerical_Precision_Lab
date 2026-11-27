# cir-numerical-precision-lab

Scientific prototype (not a production pricing library) that measures **where higher arithmetic precision
(`BigDecimal`) changes the Monte Carlo valuation of a Bermudan call on a coupon bond in the Cox-Ingersoll-Ross
model, compared with `double`** - and what it costs.

Two implementations that are identical in model, Euler scheme, time grid, random numbers (common random numbers,
the very same binary64 normals), cashflows, discounting rule and backward induction logic, differing only in the
numerical representation:

| | double | BigDecimal |
|---|---|---|
| state, cashflows, discount factors, sums, means | `double` | `BigDecimal` with explicit `MathContext` (default DECIMAL128: 34 digits, HALF_EVEN) |
| exp / log / sqrt | `java.lang.Math` | `ch.obermuhlner:big-math` 2.3.2 with guard digits, no round trip through `double` |
| random input | shared binary64 normals | exact conversion `new BigDecimal(z)` |
| regression (LSM) | Householder QR | Householder QR, same operation order |

Continuation values are estimated with **quantile clusters** of the short rate (primary method) or
**Longstaff-Schwartz** least squares (`1, x, x^2`). The double-vs-BigDecimal comparison runs in two modes:
**A** (shared clusters: identical path indices per cluster) and **B** (separate clusters: each arithmetic sorts its
own states). Policies are trained on training paths and applied out-of-sample to independent valuation paths.

Documentation:

* [docs/METHODOLOGY.md](docs/METHODOLOGY.md) - mathematics, conventions, design decisions, precision strategy, metrics, known limits
* [docs/PLAN.md](docs/PLAN.md) - development plan, phase status, acceptance criteria, risks, open decisions
* [docs/RESULTS.md](docs/RESULTS.md) - measured results of the experiments in this repository (with experiment ids)
* [CIR_Numerical_Precision_Lab.md](CIR_Numerical_Precision_Lab.md) - original specification (German)

## Build and test

Requirements: JDK 21, Maven 3.8+. Internet access for the first dependency download.

```bash
mvn -B package            # compiles, runs the 47 deterministic tests, builds target/cir-numerical-precision-lab.jar
```

## Run

```bash
# full comparison: console report + CSV under results/<experiment-id>/
java -jar target/cir-numerical-precision-lab.jar compare config/standard.properties

# quick smoke test (seconds, no CSV)
java -jar target/cir-numerical-precision-lab.jar compare config/smoke.properties

# any key can be overridden on the command line
java -jar target/cir-numerical-precision-lab.jar compare config/standard.properties decimal.precision=64 cluster.count=20

# European control cases vs closed form / Jamshidian / finite differences, Euler step series
java -jar target/cir-numerical-precision-lab.jar european config/standard.properties valuation.paths=5000

# error budget: quadrature/FD Bermudan reference, Euler step series (double), precision effect for scale
java -jar target/cir-numerical-precision-lab.jar reference config/standard.properties simulation.paths=20000 valuation.paths=20000

# accuracy of the deterministic references: quadrature settings, FD schemes, separate space/time/rMax refinement (~1.5 min)
java -jar target/cir-numerical-precision-lab.jar fdstudy config/standard.properties

# clean cluster-vs-LSM comparison: trained policies valued exactly, 40 seeds x 7 path counts x exact/Euler paths (~15 min)
java -jar target/cir-numerical-precision-lab.jar methodstudy config/standard.properties
# smaller: methodstudy.paths=500,4000 methodstudy.replications=8 methodstudy.sources=EXACT

# one-factor series (no Cartesian products); "a+b" sets several keys to the same value
java -jar target/cir-numerical-precision-lab.jar sensitivity config/standard.properties decimal.precision 16,20,34,50,64,100
java -jar target/cir-numerical-precision-lab.jar sensitivity config/standard.properties "simulation.paths+valuation.paths" 1000,4000,16000

# portfolio: independent valuation vs shared market simulation
java -jar target/cir-numerical-precision-lab.jar portfolio config/standard.properties 100

# JMH micro benchmarks (separate measurement areas, precision 16/34/64), with GC/allocation profiler
java -cp target/cir-numerical-precision-lab.jar org.openjdk.jmh.Main lab.cir.bench -prof gc
```

Profiles: `config/smoke.properties` (seconds), `config/standard.properties` (spec 19.1, ~5 s),
`config/research.properties` (20 000 + 20 000 paths, 40 clusters, 4 threads; start explicitly, minutes).
`results/run_batch1.sh` reproduces the experiment batch reported in `docs/RESULTS.md`.

## Configuration keys

| key | default | meaning |
|---|---|---|
| `cir.meanReversionSpeed`, `cir.longTermMean`, `cir.volatility`, `cir.initialRate` | - | theta, mu, sigma, r0 (decimal strings) |
| `cir.fellerMode` | `STRICT` | `STRICT` rejects `2 theta mu < sigma^2`; `ALLOW_VIOLATION` for robustness runs |
| `bond.notional`, `bond.couponRate`, `bond.maturity`, `bond.couponFrequency` | - / 1 | regular coupon bond; or `bond.cashflowTimes` + `bond.cashflowAmounts` |
| `option.strike`, `option.exerciseDates` | - | money strike; exercise dates in years, must lie on the grid |
| `simulation.stepsPerYear`, `simulation.normalsPerStep` | - / 1 | Euler grid; coupled Brownian increments (fine normals per step) |
| `simulation.paths`, `simulation.seed` | - | training paths and seed |
| `valuation.paths`, `valuation.seed` | training / seed+1000003 | independent out-of-sample paths |
| `simulation.threads` | 1 | path-parallel simulation; results are independent of the thread count |
| `cluster.count`, `cluster.minPathsPerCluster` | 10 / 20 | quantile clusters, minimum size (reduces K) |
| `continuation.methods` | `QUANTILE_CLUSTER,LEAST_SQUARES` | methods to run |
| `continuation.pathSelection` | `ALL` | `ALL` or `IN_THE_MONEY` paths in the estimate |
| `lsm.basisVariable` | `SHORT_RATE` | regression variable x in `1, x, x^2` (`SHORT_RATE` or `BOND_VALUE`) |
| `decimal.precision`, `decimal.roundingMode`, `decimal.guardDigits` | 34 / HALF_EVEN / 10 | BigDecimal working precision; guard digits for exp/log/sqrt |
| `zerobond.form` | `LOG_DIFFERENCE` | evaluation form of `ln A(tau)` (`LOG1P_RATIO` avoids a cancellation) |
| `control.kahan` | true | additional double run with Kahan-summed cluster means |
| `diagnostics.tracePaths` | 200 | paths traced along the full grid for 12.2 |
| `output.dir`, `output.csv`, `experiment.label` | results / true / - | CSV output |

## Output

Each run gets an experiment id `yyyyMMdd-HHmmss-SSS-<kind>-<label>` and a directory under `results/` with
`run_metadata.properties` (configuration, seeds, Java/JVM/CPU/OS/heap/GC, library versions, git commit) and CSV
files (`summary`, `comparisons`, `exercise_dates`, `clusters`, `propagation`, `cluster_switches`,
`regression_coefficients`, `state_differences`, `trajectory_differences`, `timings`, ...). BigDecimal values are
written as full decimal strings without conversion to double.

## Package overview

| package | responsibility |
|---|---|
| `lab.cir.model` | CIR parameters (Feller modes), cashflow schedule, Bermudan call |
| `lab.cir.analytic` | CIR zero bond (double / BigDecimal), noncentral chi-square, zero-bond and Jamshidian coupon-bond options, finite-difference and quadrature Bermudan references |
| `lab.cir.random` | AS 241 inverse normal, per-path SplitMix64 normal source (common random numbers) |
| `lab.cir.simulation` | time grid, full-truncation Euler generators (double / BigDecimal), path data |
| `lab.cir.numeric` | precision settings, `DecimalMath` (big-math facade), Householder QR (double / BigDecimal) |
| `lab.cir.induction` | quantile partition, training structure (Mode A/B), backward induction, policies, out-of-sample valuation |
| `lab.cir.stats`, `lab.cir.metrics` | sample statistics, paired differences, comparison metrics 12.1-12.11 |
| `lab.cir.experiment`, `lab.cir.report`, `lab.cir.app` | experiments, console/CSV reporting, CLI |
| `lab.cir.bench` | JMH benchmarks |

The two arithmetics are deliberately implemented as separate, line-by-line comparable classes
(`...Double` / `...Decimal`) instead of a generic numeric abstraction; only index-level logic (quantile partition,
cluster assignment, decision bookkeeping) is shared.

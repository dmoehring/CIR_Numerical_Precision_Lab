package lab.cir.induction;

import lab.cir.numeric.HouseholderQrDouble;
import lab.cir.simulation.PathDataDouble;
import lab.cir.stats.SampleStatsDouble;

import java.util.Arrays;

/**
 * Backward induction over the exercise dates in double.
 *
 * <pre>
 *   last date m:   exercise iff E_m > 0                               (V_m = E_m)
 *   j = m-1..1:    Y_j^p = exp(-(I_tau^p - I_j^p)) E_tau^p            (0 if no later exercise, tau = selected later date)
 *                  C_j   = estimate from {Y_j^p : p eligible}
 *                  exercise iff E_j^p > 0 and E_j^p > C_j^p           (tie: continue)
 *                  if exercised: tau^p := j (replaces the later cashflow)
 *   price:         mean_p exp(-I_tau^p) E_tau^p
 * </pre>
 *
 * Per path only the selected exercise index is kept; cashflows and discount factors are recomputed from
 * the stored exercise values and cumulative integrals.
 */
public final class BackwardInductionDouble {

    private BackwardInductionDouble() {
    }

    public static InductionResultDouble run(PathDataDouble d, TrainingStructure ts, InductionSettings s) {
        if (ts.source() != TrainingStructure.Source.DOUBLE) {
            throw new IllegalArgumentException("double induction requires a structure built from double states");
        }
        final int m = d.dates();
        final int paths = d.paths();
        InductionResultDouble r = new InductionResultDouble(s, ts, m, paths);
        int[] tau = new int[paths];
        Arrays.fill(tau, -1);

        // last exercise date: no continuation
        int last = m - 1;
        r.futureExerciseBefore[last] = tau.clone();
        r.target[last] = new double[paths];
        r.continuation[last] = new double[paths];
        r.exercise[last] = new boolean[paths];
        r.cluster[last] = filled(paths, -1);
        for (int p = 0; p < paths; p++) {
            if (d.exerciseValue[last][p] > 0.0) {
                r.exercise[last][p] = true;
                tau[p] = last;
            }
        }

        for (int j = m - 2; j >= 0; j--) {
            r.futureExerciseBefore[j] = tau.clone();
            double[] target = new double[paths];
            for (int p = 0; p < paths; p++) {
                int t = tau[p];
                target[p] = t < 0 ? 0.0
                        : d.exerciseValue[t][p] * Math.exp(-(d.cumulativeIntegral[t][p] - d.cumulativeIntegral[j][p]));
            }
            double[] continuation = new double[paths];
            Arrays.fill(continuation, Double.NaN);
            int[] cluster = filled(paths, -1);
            r.rules[j] = switch (s.method()) {
                case QUANTILE_CLUSTER -> fitClusters(d, ts, s, j, target, continuation, cluster);
                case LEAST_SQUARES -> fitRegression(d, ts, s, j, target, continuation);
            };
            boolean[] exercise = new boolean[paths];
            for (int p = 0; p < paths; p++) {
                double e = d.exerciseValue[j][p];
                double c = continuation[p];
                // paths outside the estimation set (NaN) are never exercised
                if (e > 0.0 && e > c) {
                    exercise[p] = true;
                    tau[p] = j;
                }
            }
            r.target[j] = target;
            r.continuation[j] = continuation;
            r.cluster[j] = cluster;
            r.exercise[j] = exercise;
        }

        for (int p = 0; p < paths; p++) {
            int t = tau[p];
            r.exerciseIndex[p] = t;
            r.pathValue[p] = t < 0 ? 0.0 : d.exerciseValue[t][p] * Math.exp(-d.cumulativeIntegral[t][p]);
        }
        r.price = SampleStatsDouble.of(r.pathValue);
        return r;
    }

    private static ClusterRuleDouble fitClusters(PathDataDouble d, TrainingStructure ts, InductionSettings s, int j,
                                                 double[] target, double[] continuation, int[] cluster) {
        QuantilePartition qp = ts.partition(j);
        int k = qp.clusterCount();
        double[] mean = new double[k];
        int[] size = new int[k];
        double[] minState = new double[k];
        double[] maxState = new double[k];
        double[] sd = new double[k];
        double[] state = d.rate[j];
        for (int c = 0; c < k; c++) {
            int start = qp.clusterStart(c);
            int end = qp.clusterEnd(c);
            int n = end - start;
            double sum = 0.0;
            double comp = 0.0; // Kahan compensation
            for (int pos = start; pos < end; pos++) {
                double y = target[qp.sortedPath(pos)];
                if (s.doubleSummation() == SummationMode.KAHAN) {
                    double yc = y - comp;
                    double t = sum + yc;
                    comp = (t - sum) - yc;
                    sum = t;
                } else {
                    sum += y;
                }
            }
            mean[c] = sum / n;
            size[c] = n;
            double ss = 0.0;
            double lo = Double.POSITIVE_INFINITY;
            double hi = Double.NEGATIVE_INFINITY;
            for (int pos = start; pos < end; pos++) {
                int p = qp.sortedPath(pos);
                double dev = target[p] - mean[c];
                ss += dev * dev;
                lo = Math.min(lo, state[p]);
                hi = Math.max(hi, state[p]);
                continuation[p] = mean[c];
                cluster[p] = c;
            }
            sd[c] = n > 1 ? Math.sqrt(ss / (n - 1)) : 0.0;
            minState[c] = lo;
            maxState[c] = hi;
        }
        return new ClusterRuleDouble(ts.upperBoundsDouble(j), mean, size, minState, maxState, sd);
    }

    private static RegressionRuleDouble fitRegression(PathDataDouble d, TrainingStructure ts, InductionSettings s, int j,
                                                      double[] target, double[] continuation) {
        int[] eligible = ts.eligible(j);
        int n = eligible.length;
        double[] beta;
        if (n < 3) {
            // too few paths for three coefficients: constant estimate (mean of targets, 0 if none)
            double sum = 0.0;
            for (int p : eligible) sum += target[p];
            beta = new double[]{n == 0 ? 0.0 : sum / n, 0.0, 0.0};
        } else {
            double[][] x = new double[3][n];
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                int p = eligible[i];
                double v = s.regressionVariable() == RegressionVariable.SHORT_RATE ? d.rate[j][p] : d.bondValue[j][p];
                x[0][i] = 1.0;
                x[1][i] = v;
                x[2][i] = v * v;
                y[i] = target[p];
            }
            beta = HouseholderQrDouble.solve(x, y);
        }
        RegressionRuleDouble rule = new RegressionRuleDouble(s.regressionVariable(), beta, n);
        for (int p : eligible) {
            continuation[p] = rule.predict(d.rate[j][p], d.bondValue[j][p]);
        }
        return rule;
    }

    static int[] filled(int n, int value) {
        int[] a = new int[n];
        Arrays.fill(a, value);
        return a;
    }
}

package lab.cir.induction;

import lab.cir.numeric.DecimalMath;
import lab.cir.numeric.HouseholderQrDecimal;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.stats.SampleStatsDecimal;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;

/**
 * Backward induction in BigDecimal; identical logic and summation order as {@link BackwardInductionDouble}.
 *
 * <p>The training structure decides the comparison mode: built from the double states it gives Mode A
 * (shared clusters / shared estimation set), built from the BigDecimal states Mode B (separate clusters).
 */
public final class BackwardInductionDecimal {

    private BackwardInductionDecimal() {
    }

    public static InductionResultDecimal run(PathDataDecimal d, TrainingStructure ts, InductionSettings s, DecimalMath math) {
        final MathContext mc = math.mc();
        final int m = d.dates();
        final int paths = d.paths();
        InductionResultDecimal r = new InductionResultDecimal(s, ts, m, paths);
        int[] tau = new int[paths];
        Arrays.fill(tau, -1);

        int last = m - 1;
        r.futureExerciseBefore[last] = tau.clone();
        r.target[last] = filled(paths, BigDecimal.ZERO);
        r.continuation[last] = filled(paths, BigDecimal.ZERO);
        r.exercise[last] = new boolean[paths];
        r.cluster[last] = BackwardInductionDouble.filled(paths, -1);
        for (int p = 0; p < paths; p++) {
            if (d.exerciseValue[last][p].signum() > 0) {
                r.exercise[last][p] = true;
                tau[p] = last;
            }
        }

        for (int j = m - 2; j >= 0; j--) {
            r.futureExerciseBefore[j] = tau.clone();
            BigDecimal[] target = new BigDecimal[paths];
            for (int p = 0; p < paths; p++) {
                int t = tau[p];
                if (t < 0) {
                    target[p] = BigDecimal.ZERO;
                } else {
                    BigDecimal exponent = d.cumulativeIntegral[t][p].subtract(d.cumulativeIntegral[j][p], mc).negate();
                    target[p] = d.exerciseValue[t][p].multiply(math.exp(exponent), mc);
                }
            }
            BigDecimal[] continuation = new BigDecimal[paths];
            int[] cluster = BackwardInductionDouble.filled(paths, -1);
            r.rules[j] = switch (s.method()) {
                case QUANTILE_CLUSTER -> fitClusters(d, ts, j, target, continuation, cluster, math);
                case LEAST_SQUARES -> fitRegression(d, ts, s, j, target, continuation, math);
            };
            boolean[] exercise = new boolean[paths];
            for (int p = 0; p < paths; p++) {
                BigDecimal e = d.exerciseValue[j][p];
                BigDecimal c = continuation[p];
                if (c != null && e.signum() > 0 && e.compareTo(c) > 0) {
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
            r.pathValue[p] = t < 0 ? BigDecimal.ZERO
                    : d.exerciseValue[t][p].multiply(math.exp(d.cumulativeIntegral[t][p].negate()), mc);
        }
        r.price = SampleStatsDecimal.of(r.pathValue, math);
        return r;
    }

    private static ClusterRuleDecimal fitClusters(PathDataDecimal d, TrainingStructure ts, int j, BigDecimal[] target,
                                                  BigDecimal[] continuation, int[] cluster, DecimalMath math) {
        MathContext mc = math.mc();
        QuantilePartition qp = ts.partition(j);
        int k = qp.clusterCount();
        BigDecimal[] mean = new BigDecimal[k];
        int[] size = new int[k];
        BigDecimal[] minState = new BigDecimal[k];
        BigDecimal[] maxState = new BigDecimal[k];
        BigDecimal[] sd = new BigDecimal[k];
        BigDecimal[] state = d.rate[j];
        for (int c = 0; c < k; c++) {
            int start = qp.clusterStart(c);
            int end = qp.clusterEnd(c);
            int n = end - start;
            BigDecimal sum = BigDecimal.ZERO;
            for (int pos = start; pos < end; pos++) {
                sum = sum.add(target[qp.sortedPath(pos)], mc);
            }
            mean[c] = sum.divide(BigDecimal.valueOf(n), mc);
            size[c] = n;
            BigDecimal ss = BigDecimal.ZERO;
            BigDecimal lo = null;
            BigDecimal hi = null;
            for (int pos = start; pos < end; pos++) {
                int p = qp.sortedPath(pos);
                BigDecimal dev = target[p].subtract(mean[c], mc);
                ss = ss.add(dev.multiply(dev, mc), mc);
                if (lo == null || state[p].compareTo(lo) < 0) lo = state[p];
                if (hi == null || state[p].compareTo(hi) > 0) hi = state[p];
                continuation[p] = mean[c];
                cluster[p] = c;
            }
            sd[c] = n > 1 ? math.sqrt(ss.divide(BigDecimal.valueOf(n - 1L), mc)) : BigDecimal.ZERO;
            minState[c] = lo;
            maxState[c] = hi;
        }
        boolean shared = ts.source() == TrainingStructure.Source.DOUBLE;
        return new ClusterRuleDecimal(shared ? ts.upperBoundsDouble(j) : null,
                shared ? null : ts.upperBoundsDecimal(j), mean, size, minState, maxState, sd);
    }

    private static RegressionRuleDecimal fitRegression(PathDataDecimal d, TrainingStructure ts, InductionSettings s, int j,
                                                       BigDecimal[] target, BigDecimal[] continuation, DecimalMath math) {
        MathContext mc = math.mc();
        int[] eligible = ts.eligible(j);
        int n = eligible.length;
        BigDecimal[] beta;
        if (n < 3) {
            BigDecimal sum = BigDecimal.ZERO;
            for (int p : eligible) sum = sum.add(target[p], mc);
            beta = new BigDecimal[]{n == 0 ? BigDecimal.ZERO : sum.divide(BigDecimal.valueOf(n), mc),
                    BigDecimal.ZERO, BigDecimal.ZERO};
        } else {
            BigDecimal[][] x = new BigDecimal[3][n];
            BigDecimal[] y = new BigDecimal[n];
            for (int i = 0; i < n; i++) {
                int p = eligible[i];
                BigDecimal v = s.regressionVariable() == RegressionVariable.SHORT_RATE ? d.rate[j][p] : d.bondValue[j][p];
                x[0][i] = BigDecimal.ONE;
                x[1][i] = v;
                x[2][i] = v.multiply(v, mc);
                y[i] = target[p];
            }
            beta = HouseholderQrDecimal.solve(x, y, math);
        }
        RegressionRuleDecimal rule = new RegressionRuleDecimal(s.regressionVariable(), beta, n);
        for (int p : eligible) {
            continuation[p] = rule.predict(d.rate[j][p], d.bondValue[j][p], mc);
        }
        return rule;
    }

    private static BigDecimal[] filled(int n, BigDecimal v) {
        BigDecimal[] a = new BigDecimal[n];
        Arrays.fill(a, v);
        return a;
    }
}

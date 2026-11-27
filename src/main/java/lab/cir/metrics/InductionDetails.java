package lab.cir.metrics;

import lab.cir.induction.ClusterRuleDecimal;
import lab.cir.induction.ClusterRuleDouble;
import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.QuantilePartition;
import lab.cir.induction.RegressionRuleDecimal;
import lab.cir.induction.RegressionRuleDouble;
import lab.cir.numeric.DecimalMath;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.stats.DifferenceStats;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;

/**
 * In-sample details of a double vs BigDecimal backward induction: cluster continuation values (12.4),
 * cluster switches near quantile boundaries (12.5), propagation through the backward induction (12.8),
 * regression coefficients (LSM).
 */
public final class InductionDetails {

    /**
     * One cluster (index k) at one date. {@code standardErrorOfMean = sd/sqrt(n)} of the double targets is the
     * statistical yardstick for {@code |C_dec - C_dbl|}.
     */
    public record ClusterRow(int date, int cluster, int sizeDouble, int sizeDecimal, double meanDouble, BigDecimal meanDecimal,
                             BigDecimal absDifference, double relDifference, double minState, double maxState,
                             double targetStdDev, double standardErrorOfMean, double differenceOverStandardError) {
    }

    /** Propagation (12.8): later decision differences change the targets of earlier dates. */
    public record PropagationRow(int date, int differentFutureExercise, DifferenceStats targetDiff,
                                 int laterDecisionDifferences) {
    }

    /** Cluster switches in Mode B (12.5); rank distances refer to the double partition. */
    public record SwitchRow(int date, int same, int switched, int switchedWithinOneRank, int maxRankDistance,
                            double maxRateDifference, double meanRateDifference) {
    }

    public record CoefficientRow(int date, int coefficient, double valueDouble, BigDecimal valueDecimal, double relDifference) {
    }

    public final List<ClusterRow> clusters = new ArrayList<>();
    public final List<PropagationRow> propagation = new ArrayList<>();
    public final List<SwitchRow> switches = new ArrayList<>();
    public final List<CoefficientRow> coefficients = new ArrayList<>();

    public InductionDetails(InductionResultDouble a, InductionResultDecimal b, PathDataDouble dd, PathDataDecimal dc,
                            DecimalMath math) {
        MathContext mc = math.mc();
        int m = a.dates();
        int paths = a.paths();
        for (int j = 0; j < m - 1; j++) {
            // 12.4 cluster table / coefficients
            if (a.rules[j] instanceof ClusterRuleDouble cd && b.rules[j] instanceof ClusterRuleDecimal cb) {
                int k = Math.max(cd.clusterCount(), cb.clusterCount());
                for (int c = 0; c < k; c++) {
                    boolean hasD = c < cd.clusterCount();
                    boolean hasB = c < cb.clusterCount();
                    double meanD = hasD ? cd.mean()[c] : Double.NaN;
                    BigDecimal meanB = hasB ? cb.mean()[c] : null;
                    BigDecimal diff = hasD && hasB ? meanB.subtract(new BigDecimal(meanD), mc).abs() : null;
                    double rel = diff == null || meanD == 0.0 ? Double.NaN : diff.doubleValue() / Math.abs(meanD);
                    double sd = hasD ? cd.targetStdDev()[c] : Double.NaN;
                    double se = hasD ? sd / Math.sqrt(cd.size()[c]) : Double.NaN;
                    clusters.add(new ClusterRow(j, c, hasD ? cd.size()[c] : 0, hasB ? cb.size()[c] : 0, meanD, meanB, diff, rel,
                            hasD ? cd.minState()[c] : Double.NaN, hasD ? cd.maxState()[c] : Double.NaN, sd, se,
                            diff == null || se == 0.0 ? Double.NaN : diff.doubleValue() / se));
                }
            } else if (a.rules[j] instanceof RegressionRuleDouble rd && b.rules[j] instanceof RegressionRuleDecimal rb) {
                for (int c = 0; c < rd.coefficients().length; c++) {
                    double v = rd.coefficients()[c];
                    BigDecimal w = rb.coefficients()[c];
                    double rel = v == 0.0 ? Double.NaN : w.subtract(new BigDecimal(v), mc).abs().doubleValue() / Math.abs(v);
                    coefficients.add(new CoefficientRow(j, c, v, w, rel));
                }
            }
            // 12.8 propagation
            int diffFuture = 0;
            DifferenceStats td = new DifferenceStats(1e-12, mc);
            for (int p = 0; p < paths; p++) {
                if (a.futureExerciseBefore[j][p] != b.futureExerciseBefore[j][p]) diffFuture++;
                td.add(b.target[j][p], a.target[j][p]);
            }
            int later = 0;
            for (int l = j + 1; l < m; l++) {
                for (int p = 0; p < paths; p++) {
                    if (a.exercise[l][p] != b.exercise[l][p]) later++;
                }
            }
            propagation.add(new PropagationRow(j, diffFuture, td, later));
            // 12.5 cluster switches (only meaningful if b has its own clusters)
            QuantilePartition qd = a.structure.partition(j);
            QuantilePartition qb = b.structure.partition(j);
            if (qd != null && qb != null) {
                int same = 0, switched = 0, near = 0, maxRank = 0;
                double maxDr = 0, sumDr = 0;
                for (int p = 0; p < paths; p++) {
                    int kd = qd.clusterOf(p);
                    int kb = qb.clusterOf(p);
                    if (kd < 0 || kb < 0) continue;
                    if (kd == kb) {
                        same++;
                    } else {
                        switched++;
                        int rd = qd.rankDistanceToBoundary(p);
                        if (rd <= 1) near++;
                        maxRank = Math.max(maxRank, rd);
                        double dr = dc.rate[j][p].subtract(new BigDecimal(dd.rate[j][p]), mc).abs().doubleValue();
                        maxDr = Math.max(maxDr, dr);
                        sumDr += dr;
                    }
                }
                switches.add(new SwitchRow(j, same, switched, near, switched == 0 ? 0 : maxRank, maxDr,
                        switched == 0 ? 0.0 : sumDr / switched));
            }
        }
    }
}

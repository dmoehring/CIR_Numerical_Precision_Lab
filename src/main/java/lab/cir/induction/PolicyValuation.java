package lab.cir.induction;

import lab.cir.numeric.DecimalMath;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.stats.SampleStatsDecimal;
import lab.cir.stats.SampleStatsDouble;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;

/**
 * Applies fixed trained exercise rules to independent valuation paths in chronological order
 * (out-of-sample valuation). The realised exercise is the first date with {@code E > 0 && E > C(state)};
 * at the last date {@code C = 0}.
 *
 * <p>Mode A (shared clusters): the BigDecimal policy assigns valuation paths to clusters with the double
 * boundaries applied to the double valuation states (same assignment as the double policy); only the
 * continuation values per cluster differ. Mode B: each arithmetic uses its own boundaries and states.
 */
public final class PolicyValuation {

    private PolicyValuation() {
    }

    public static ValuationResultDouble valueDouble(PathDataDouble v, InductionResultDouble trained) {
        int m = v.dates();
        int paths = v.paths();
        ValuationResultDouble r = new ValuationResultDouble(m, paths);
        for (int j = 0; j < m; j++) {
            Arrays.fill(r.cluster[j], -1);
            DateRuleDouble rule = trained.rules[j];
            for (int p = 0; p < paths; p++) {
                double c;
                if (rule == null) {
                    c = 0.0;
                } else if (rule instanceof ClusterRuleDouble cr) {
                    if (cr.clusterCount() == 0) {
                        c = Double.NaN;
                    } else {
                        int k = cr.assign(v.rate[j][p]);
                        r.cluster[j][p] = k;
                        c = cr.mean()[k];
                    }
                } else {
                    c = ((RegressionRuleDouble) rule).predict(v.rate[j][p], v.bondValue[j][p]);
                }
                r.continuation[j][p] = c;
                double e = v.exerciseValue[j][p];
                r.decision[j][p] = e > 0.0 && e > c;
            }
        }
        for (int p = 0; p < paths; p++) {
            int t = firstExercise(r.decision, p);
            r.exerciseIndex[p] = t;
            r.pathValue[p] = t < 0 ? 0.0 : v.exerciseValue[t][p] * Math.exp(-v.cumulativeIntegral[t][p]);
        }
        r.price = SampleStatsDouble.of(r.pathValue);
        return r;
    }

    /**
     * @param sharedDoubleStates double valuation paths (same random numbers); required when the trained policy uses
     *                           shared cluster assignment (Mode A), ignored otherwise
     */
    public static ValuationResultDecimal valueDecimal(PathDataDecimal v, InductionResultDecimal trained,
                                                      PathDataDouble sharedDoubleStates, DecimalMath math) {
        MathContext mc = math.mc();
        int m = v.dates();
        int paths = v.paths();
        ValuationResultDecimal r = new ValuationResultDecimal(m, paths);
        for (int j = 0; j < m; j++) {
            Arrays.fill(r.cluster[j], -1);
            DateRuleDecimal rule = trained.rules[j];
            for (int p = 0; p < paths; p++) {
                BigDecimal c;
                if (rule == null) {
                    c = BigDecimal.ZERO;
                } else if (rule instanceof ClusterRuleDecimal cr) {
                    if (cr.clusterCount() == 0) {
                        c = null;
                    } else {
                        int k = cr.sharedAssignment() ? cr.assignShared(sharedDoubleStates.rate[j][p]) : cr.assign(v.rate[j][p]);
                        r.cluster[j][p] = k;
                        c = cr.mean()[k];
                    }
                } else {
                    c = ((RegressionRuleDecimal) rule).predict(v.rate[j][p], v.bondValue[j][p], mc);
                }
                r.continuation[j][p] = c;
                BigDecimal e = v.exerciseValue[j][p];
                r.decision[j][p] = c != null && e.signum() > 0 && e.compareTo(c) > 0;
            }
        }
        for (int p = 0; p < paths; p++) {
            int t = firstExercise(r.decision, p);
            r.exerciseIndex[p] = t;
            r.pathValue[p] = t < 0 ? BigDecimal.ZERO
                    : v.exerciseValue[t][p].multiply(math.exp(v.cumulativeIntegral[t][p].negate()), mc);
        }
        r.price = SampleStatsDecimal.of(r.pathValue, math);
        return r;
    }

    /**
     * Frozen-policy control: BigDecimal cashflows for externally fixed exercise indices (e.g. those of the double
     * policy). Separates pure cashflow arithmetic from the effect of differing decisions.
     */
    public static BigDecimal[] decimalCashflowsForFixedExercise(PathDataDecimal v, int[] exerciseIndex, DecimalMath math) {
        MathContext mc = math.mc();
        BigDecimal[] values = new BigDecimal[v.paths()];
        for (int p = 0; p < values.length; p++) {
            int t = exerciseIndex[p];
            values[p] = t < 0 ? BigDecimal.ZERO
                    : v.exerciseValue[t][p].multiply(math.exp(v.cumulativeIntegral[t][p].negate()), mc);
        }
        return values;
    }

    private static int firstExercise(boolean[][] decision, int p) {
        for (int j = 0; j < decision.length; j++) {
            if (decision[j][p]) return j;
        }
        return -1;
    }
}

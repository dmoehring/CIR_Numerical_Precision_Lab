package lab.cir.metrics;

import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.ValuationResultDecimal;
import lab.cir.induction.ValuationResultDouble;
import lab.cir.stats.SampleStatsDecimal;
import lab.cir.stats.SampleStatsDouble;

import java.math.BigDecimal;

/** Uniform read-only view on in-sample (backward induction) and out-of-sample (forward valuation) results. */
public final class PolicyView {

    private PolicyView() {
    }

    public record OfDouble(String kind, boolean[][] decision, double[][] continuation, int[][] cluster,
                           int[] exerciseIndex, double[] pathValue, SampleStatsDouble price) {
        public static OfDouble of(InductionResultDouble r) {
            return new OfDouble("in-sample", r.exercise, r.continuation, r.cluster, r.exerciseIndex, r.pathValue, r.price);
        }

        public static OfDouble of(ValuationResultDouble r) {
            return new OfDouble("out-of-sample", r.decision, r.continuation, r.cluster, r.exerciseIndex, r.pathValue, r.price);
        }
    }

    public record OfDecimal(String kind, boolean[][] decision, BigDecimal[][] continuation, int[][] cluster,
                            int[] exerciseIndex, BigDecimal[] pathValue, SampleStatsDecimal price) {
        public static OfDecimal of(InductionResultDecimal r) {
            return new OfDecimal("in-sample", r.exercise, r.continuation, r.cluster, r.exerciseIndex, r.pathValue, r.price);
        }

        public static OfDecimal of(ValuationResultDecimal r) {
            return new OfDecimal("out-of-sample", r.decision, r.continuation, r.cluster, r.exerciseIndex, r.pathValue, r.price);
        }
    }
}

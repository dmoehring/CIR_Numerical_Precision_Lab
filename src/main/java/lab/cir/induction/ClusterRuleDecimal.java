package lab.cir.induction;

import java.math.BigDecimal;

/**
 * BigDecimal quantile cluster rule. In Mode A (shared clusters) the assignment of new states is done with the
 * <b>double</b> boundaries applied to the double states ({@code sharedUpperBounds != null}); in Mode B with the
 * own BigDecimal boundaries.
 */
public record ClusterRuleDecimal(double[] sharedUpperBounds, BigDecimal[] upperBounds, BigDecimal[] mean, int[] size,
                                 BigDecimal[] minState, BigDecimal[] maxState, BigDecimal[] targetStdDev)
        implements DateRuleDecimal {

    public int clusterCount() { return mean.length; }

    public boolean sharedAssignment() { return sharedUpperBounds != null; }

    public int assign(BigDecimal state) {
        int lo = 0;
        int hi = upperBounds.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (state.compareTo(upperBounds[mid]) <= 0) hi = mid; else lo = mid + 1;
        }
        return lo;
    }

    public int assignShared(double doubleState) {
        return ClusterRuleDouble.assign(sharedUpperBounds, doubleState);
    }
}

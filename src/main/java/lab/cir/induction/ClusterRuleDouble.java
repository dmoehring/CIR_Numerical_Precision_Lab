package lab.cir.induction;

/**
 * Quantile cluster rule: cluster {@code k} covers states {@code (upperBounds[k-1], upperBounds[k]]},
 * the first cluster extends to -infinity and the last to +infinity. {@code mean[k]} is the estimated
 * continuation value; the other arrays are diagnostics of the training clusters.
 */
public record ClusterRuleDouble(double[] upperBounds, double[] mean, int[] size, double[] minState,
                                double[] maxState, double[] targetStdDev) implements DateRuleDouble {

    public int clusterCount() { return mean.length; }

    /** Cluster of a new state: smallest k with {@code state <= upperBounds[k]}, else the last cluster. */
    public int assign(double state) {
        return assign(upperBounds, state);
    }

    public static int assign(double[] upperBounds, double state) {
        int lo = 0;
        int hi = upperBounds.length; // answer in [0, K-1]
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (state <= upperBounds[mid]) hi = mid; else lo = mid + 1;
        }
        return lo;
    }
}

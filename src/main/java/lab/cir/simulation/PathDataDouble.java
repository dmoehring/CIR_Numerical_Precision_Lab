package lab.cir.simulation;

/**
 * Simulated states at the exercise dates, layout {@code [exerciseIndex][pathIndex]}.
 * Only what the valuation needs is kept (no full paths).
 */
public final class PathDataDouble {

    /** Effective (non-negative) short rate {@code x^+} at {@code e_j}; used for bond values, clustering, regression. */
    public final double[][] rate;
    /** Left-rectangle integral {@code sum_{k < k_j} x_k^+ dt} of the effective rate from 0 to {@code e_j}. */
    public final double[][] cumulativeIntegral;
    /** Remaining bond value {@code B(e_j)}. */
    public final double[][] bondValue;
    /** Immediate exercise value {@code (B(e_j) - K)^+}. */
    public final double[][] exerciseValue;
    /** Number of grid points (k >= 1) at which the internal Euler state was negative, per path. */
    public final int[] negativeStateCount;

    public PathDataDouble(int dates, int paths) {
        rate = new double[dates][paths];
        cumulativeIntegral = new double[dates][paths];
        bondValue = new double[dates][paths];
        exerciseValue = new double[dates][paths];
        negativeStateCount = new int[paths];
    }

    /** View sharing rate and integral rows of {@code source} (selected observation dates) with fresh value rows. */
    public PathDataDouble(PathDataDouble source, int[] observationIndices) {
        int m = observationIndices.length;
        int paths = source.paths();
        rate = new double[m][];
        cumulativeIntegral = new double[m][];
        bondValue = new double[m][paths];
        exerciseValue = new double[m][paths];
        for (int j = 0; j < m; j++) {
            rate[j] = source.rate[observationIndices[j]];
            cumulativeIntegral[j] = source.cumulativeIntegral[observationIndices[j]];
        }
        negativeStateCount = source.negativeStateCount;
    }

    public int paths() { return negativeStateCount.length; }
    public int dates() { return rate.length; }

    public long totalNegativeStates() {
        long s = 0;
        for (int c : negativeStateCount) s += c;
        return s;
    }

    public int pathsWithNegativeStates() {
        int s = 0;
        for (int c : negativeStateCount) if (c > 0) s++;
        return s;
    }
}

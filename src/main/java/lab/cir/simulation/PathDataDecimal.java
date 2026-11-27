package lab.cir.simulation;

import java.math.BigDecimal;

/** BigDecimal counterpart of {@link PathDataDouble}. */
public final class PathDataDecimal {

    public final BigDecimal[][] rate;
    public final BigDecimal[][] cumulativeIntegral;
    public final BigDecimal[][] bondValue;
    public final BigDecimal[][] exerciseValue;
    public final int[] negativeStateCount;

    public PathDataDecimal(int dates, int paths) {
        rate = new BigDecimal[dates][paths];
        cumulativeIntegral = new BigDecimal[dates][paths];
        bondValue = new BigDecimal[dates][paths];
        exerciseValue = new BigDecimal[dates][paths];
        negativeStateCount = new int[paths];
    }

    /** View sharing rate and integral rows of {@code source} (selected observation dates) with fresh value rows. */
    public PathDataDecimal(PathDataDecimal source, int[] observationIndices) {
        int m = observationIndices.length;
        int paths = source.paths();
        rate = new BigDecimal[m][];
        cumulativeIntegral = new BigDecimal[m][];
        bondValue = new BigDecimal[m][paths];
        exerciseValue = new BigDecimal[m][paths];
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

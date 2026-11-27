package lab.cir.stats;

/**
 * Sample mean, standard deviation (n-1), standard error and normal 95% confidence interval in double.
 * Sums are naive and in index order (two-pass variance).
 */
public record SampleStatsDouble(int n, double mean, double stdDev, double standardError) {

    public static final double Z_95 = 1.959963984540054;

    public static SampleStatsDouble of(double[] values) {
        int n = values.length;
        if (n == 0) {
            return new SampleStatsDouble(0, Double.NaN, Double.NaN, Double.NaN);
        }
        double sum = 0.0;
        for (double v : values) sum += v;
        double mean = sum / n;
        if (n == 1) {
            return new SampleStatsDouble(1, mean, Double.NaN, Double.NaN);
        }
        double ss = 0.0;
        for (double v : values) {
            double d = v - mean;
            ss += d * d;
        }
        double sd = Math.sqrt(ss / (n - 1));
        return new SampleStatsDouble(n, mean, sd, sd / Math.sqrt(n));
    }

    public double ciLow() { return mean - Z_95 * standardError; }
    public double ciHigh() { return mean + Z_95 * standardError; }
}

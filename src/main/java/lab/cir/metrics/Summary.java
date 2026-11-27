package lab.cir.metrics;

import java.util.Arrays;

/** Median / mean / max of a list of non-negative diagnostic values. */
public record Summary(int count, double median, double mean, double max) {

    public static final Summary EMPTY = new Summary(0, Double.NaN, Double.NaN, Double.NaN);

    public static Summary of(double[] values, int count) {
        double[] v = Arrays.stream(values, 0, count).filter(x -> !Double.isNaN(x)).sorted().toArray();
        count = v.length;
        if (count == 0) {
            return EMPTY;
        }
        double sum = 0;
        for (double x : v) sum += x;
        double median = count % 2 == 1 ? v[count / 2] : 0.5 * (v[count / 2 - 1] + v[count / 2]);
        return new Summary(count, median, sum / count, v[count - 1]);
    }
}

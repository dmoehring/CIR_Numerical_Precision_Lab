package lab.cir.stats;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Element-wise differences {@code d_i = decimal_i - double_i}, computed exactly in BigDecimal
 * ({@code new BigDecimal(double)} is exact) and rounded to the given MathContext.
 *
 * <p>Reports: count, mean |d|, max |d|, mean d, RMSE, maximum relative difference |d| / |decimal|
 * (only over references with |decimal| > relativeFloor; the number of references at or below the floor is
 * reported separately), and the distance in units in the last place of the double value.
 */
public final class DifferenceStats {

    private int n;
    private int nonZero;
    private int nearZeroReferences;
    private double sumAbs;
    private double sumSigned;
    private double sumSquares;
    private double maxAbs;
    private double maxRel;
    private double sumUlp;
    private double maxUlp;
    private final double relativeFloor;
    private final MathContext mc;

    public DifferenceStats(double relativeFloor, MathContext mc) {
        this.relativeFloor = relativeFloor;
        this.mc = mc;
    }

    /** Accumulates one pair; returns the signed difference decimal - double. */
    public BigDecimal add(BigDecimal decimal, double dbl) {
        BigDecimal diff = decimal.subtract(new BigDecimal(dbl), mc);
        // the aggregated statistics are diagnostics only; |d| is far below double range limits
        double d = diff.doubleValue();
        double a = Math.abs(d);
        n++;
        if (diff.signum() != 0) nonZero++;
        sumAbs += a;
        sumSigned += d;
        sumSquares += d * d;
        maxAbs = Math.max(maxAbs, a);
        double ref = Math.abs(decimal.doubleValue());
        if (ref > relativeFloor) {
            maxRel = Math.max(maxRel, a / ref);
        } else {
            nearZeroReferences++;
        }
        double ulp = Math.ulp(dbl);
        double u = a / ulp;
        sumUlp += u;
        maxUlp = Math.max(maxUlp, u);
        return diff;
    }

    public int count() { return n; }
    public int nonZeroCount() { return nonZero; }
    public int nearZeroReferenceCount() { return nearZeroReferences; }
    public double meanAbs() { return n == 0 ? Double.NaN : sumAbs / n; }
    public double meanSigned() { return n == 0 ? Double.NaN : sumSigned / n; }
    public double rmse() { return n == 0 ? Double.NaN : Math.sqrt(sumSquares / n); }
    public double maxAbs() { return maxAbs; }
    public double maxRel() { return maxRel; }
    public double meanUlp() { return n == 0 ? Double.NaN : sumUlp / n; }
    public double maxUlp() { return maxUlp; }
}

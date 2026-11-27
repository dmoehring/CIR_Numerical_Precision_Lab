package lab.cir.stats;

import lab.cir.numeric.DecimalMath;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Paired analysis of path-wise discounted option cashflows {@code d_p = Y0_p^dec - Y0_p^dbl}
 * (common random numbers): mean, standard error {@code s_d / sqrt(n)} and normal 95% CI of the mean difference,
 * plus max |d|, RMSE and the number of paths with a non-zero difference. Computed in BigDecimal.
 */
public record PairedDifference(int n, BigDecimal meanDifference, BigDecimal standardError,
                               BigDecimal maxAbsDifference, BigDecimal rmse, int nonZeroPaths,
                               int pathsWithDifferentExercise) {

    public static PairedDifference of(BigDecimal[] decimal, double[] dbl, int[] exerciseDecimal, int[] exerciseDouble,
                                      DecimalMath math) {
        MathContext mc = math.mc();
        int n = dbl.length;
        BigDecimal[] d = new BigDecimal[n];
        BigDecimal maxAbs = BigDecimal.ZERO;
        BigDecimal sumSq = BigDecimal.ZERO;
        int nonZero = 0;
        int differentExercise = 0;
        for (int p = 0; p < n; p++) {
            d[p] = decimal[p].subtract(new BigDecimal(dbl[p]), mc);
            if (d[p].signum() != 0) nonZero++;
            BigDecimal a = d[p].abs();
            if (a.compareTo(maxAbs) > 0) maxAbs = a;
            sumSq = sumSq.add(d[p].multiply(d[p], mc), mc);
            if (exerciseDecimal != null && exerciseDecimal[p] != exerciseDouble[p]) differentExercise++;
        }
        SampleStatsDecimal s = SampleStatsDecimal.of(d, math);
        BigDecimal rmse = math.sqrt(sumSq.divide(BigDecimal.valueOf(n), mc));
        return new PairedDifference(n, s.mean(), s.standardError(), maxAbs, rmse, nonZero, differentExercise);
    }

    public BigDecimal ciLow(MathContext mc) {
        return meanDifference.subtract(new BigDecimal(SampleStatsDouble.Z_95).multiply(standardError, mc), mc);
    }

    public BigDecimal ciHigh(MathContext mc) {
        return meanDifference.add(new BigDecimal(SampleStatsDouble.Z_95).multiply(standardError, mc), mc);
    }
}

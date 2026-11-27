package lab.cir.stats;

import lab.cir.numeric.DecimalMath;

import java.math.BigDecimal;
import java.math.MathContext;

/** BigDecimal counterpart of {@link SampleStatsDouble} (all arithmetic at the working MathContext). */
public record SampleStatsDecimal(int n, BigDecimal mean, BigDecimal stdDev, BigDecimal standardError) {

    public static SampleStatsDecimal of(BigDecimal[] values, DecimalMath math) {
        MathContext mc = math.mc();
        int n = values.length;
        if (n == 0) {
            return new SampleStatsDecimal(0, null, null, null);
        }
        BigDecimal sum = BigDecimal.ZERO;
        for (BigDecimal v : values) sum = sum.add(v, mc);
        BigDecimal bn = BigDecimal.valueOf(n);
        BigDecimal mean = sum.divide(bn, mc);
        if (n == 1) {
            return new SampleStatsDecimal(1, mean, null, null);
        }
        BigDecimal ss = BigDecimal.ZERO;
        for (BigDecimal v : values) {
            BigDecimal d = v.subtract(mean, mc);
            ss = ss.add(d.multiply(d, mc), mc);
        }
        BigDecimal sd = math.sqrt(ss.divide(BigDecimal.valueOf(n - 1L), mc));
        return new SampleStatsDecimal(n, mean, sd, sd.divide(math.sqrt(bn), mc));
    }

    public BigDecimal ciLow(MathContext mc) {
        return mean.subtract(new BigDecimal(SampleStatsDouble.Z_95).multiply(standardError, mc), mc);
    }

    public BigDecimal ciHigh(MathContext mc) {
        return mean.add(new BigDecimal(SampleStatsDouble.Z_95).multiply(standardError, mc), mc);
    }
}

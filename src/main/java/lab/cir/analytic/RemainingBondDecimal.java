package lab.cir.analytic;

import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;
import lab.cir.numeric.DecimalMath;

import java.math.BigDecimal;
import java.math.MathContext;

/** BigDecimal counterpart of {@link RemainingBondDouble} (same summation order). */
public final class RemainingBondDecimal {

    private final BigDecimal[][] amounts;
    private final BigDecimal[][] logA;
    private final BigDecimal[][] b;
    private final BigDecimal strike;
    private final DecimalMath math;
    private final MathContext mc;

    public RemainingBondDecimal(BermudanCallOption option, CirZeroBondDecimal zeroBond) {
        this.math = zeroBond.math();
        this.mc = math.mc();
        CashflowSchedule cf = option.underlying();
        int m = option.exerciseDateCount();
        this.amounts = new BigDecimal[m][];
        this.logA = new BigDecimal[m][];
        this.b = new BigDecimal[m][];
        this.strike = option.strike();
        for (int j = 0; j < m; j++) {
            int first = cf.firstIndexAfter(option.exerciseDate(j));
            int n = cf.size() - first;
            amounts[j] = new BigDecimal[n];
            logA[j] = new BigDecimal[n];
            b[j] = new BigDecimal[n];
            BigDecimal e = option.exerciseDate(j);
            for (int i = 0; i < n; i++) {
                amounts[j][i] = cf.amount(first + i);
                // time to payment is exact in decimal (inputs are decimal strings)
                CirZeroBondDecimal.Coefficients c = zeroBond.coefficients(cf.time(first + i).subtract(e));
                logA[j][i] = c.logA();
                b[j][i] = c.b();
            }
        }
    }

    public BigDecimal bondValue(int exerciseIndex, BigDecimal rate) {
        BigDecimal[] a = amounts[exerciseIndex];
        BigDecimal[] la = logA[exerciseIndex];
        BigDecimal[] bb = b[exerciseIndex];
        BigDecimal sum = BigDecimal.ZERO;
        for (int i = 0; i < a.length; i++) {
            BigDecimal df = math.exp(la[i].subtract(bb[i].multiply(rate, mc), mc));
            sum = sum.add(a[i].multiply(df, mc), mc);
        }
        return sum;
    }

    public BigDecimal exerciseValue(BigDecimal bondValue) {
        BigDecimal v = bondValue.subtract(strike, mc);
        return v.signum() > 0 ? v : BigDecimal.ZERO;
    }

    public BigDecimal strike() { return strike; }
}

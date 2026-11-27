package lab.cir.analytic;

import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;

/**
 * Value of the remaining cashflow stream at each exercise date,
 * {@code B(e_j, r) = sum_{t_i > e_j} C_i P(e_j, t_i; r)}, with precomputed zero-bond coefficients
 * (they depend only on {@code t_i - e_j}, not on the simulated state). Summation in increasing {@code i}.
 */
public final class RemainingBondDouble {

    private final double[][] amounts;
    private final double[][] logA;
    private final double[][] b;
    private final double strike;

    public RemainingBondDouble(BermudanCallOption option, CirZeroBondDouble zeroBond) {
        CashflowSchedule cf = option.underlying();
        int m = option.exerciseDateCount();
        this.amounts = new double[m][];
        this.logA = new double[m][];
        this.b = new double[m][];
        this.strike = option.strikeDouble();
        for (int j = 0; j < m; j++) {
            int first = cf.firstIndexAfter(option.exerciseDate(j));
            int n = cf.size() - first;
            amounts[j] = new double[n];
            logA[j] = new double[n];
            b[j] = new double[n];
            double e = option.exerciseDateDouble(j);
            for (int i = 0; i < n; i++) {
                amounts[j][i] = cf.amountDouble(first + i);
                CirZeroBondDouble.Coefficients c = zeroBond.coefficients(cf.timeDouble(first + i) - e);
                logA[j][i] = c.logA();
                b[j][i] = c.b();
            }
        }
    }

    public double bondValue(int exerciseIndex, double rate) {
        double[] a = amounts[exerciseIndex];
        double[] la = logA[exerciseIndex];
        double[] bb = b[exerciseIndex];
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * Math.exp(la[i] - bb[i] * rate);
        }
        return sum;
    }

    /** Immediate exercise value {@code (B - K)^+}. */
    public double exerciseValue(double bondValue) {
        double v = bondValue - strike;
        return v > 0.0 ? v : 0.0;
    }

    public double strike() { return strike; }
}

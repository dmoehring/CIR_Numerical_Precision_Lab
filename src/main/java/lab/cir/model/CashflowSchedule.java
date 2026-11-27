package lab.cir.model;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * A stream of strictly positive cashflows {@code C_1..C_n} at strictly increasing times
 * {@code 0 < t_1 < ... < t_n} (in years, no day count conventions).
 *
 * <p>Amounts are money amounts in the same unit as the option strike; no implicit notional scaling.
 */
public final class CashflowSchedule {

    private final BigDecimal[] times;
    private final BigDecimal[] amounts;
    private final double[] timesDouble;
    private final double[] amountsDouble;

    public CashflowSchedule(List<BigDecimal> times, List<BigDecimal> amounts) {
        Objects.requireNonNull(times, "times");
        Objects.requireNonNull(amounts, "amounts");
        if (times.isEmpty() || times.size() != amounts.size()) {
            throw new IllegalArgumentException("times and amounts must be non-empty and of equal length");
        }
        this.times = times.toArray(BigDecimal[]::new);
        this.amounts = amounts.toArray(BigDecimal[]::new);
        for (int i = 0; i < this.times.length; i++) {
            if (this.times[i].signum() <= 0) {
                throw new IllegalArgumentException("cashflow times must be > 0: " + this.times[i]);
            }
            if (i > 0 && this.times[i].compareTo(this.times[i - 1]) <= 0) {
                throw new IllegalArgumentException("cashflow times must be strictly increasing");
            }
            if (this.amounts[i].signum() <= 0) {
                throw new IllegalArgumentException("cashflow amounts must be > 0: " + this.amounts[i]);
            }
        }
        this.timesDouble = Arrays.stream(this.times).mapToDouble(CirParameters::toDouble).toArray();
        this.amountsDouble = Arrays.stream(this.amounts).mapToDouble(CirParameters::toDouble).toArray();
    }

    /**
     * Regular coupon bond: coupons {@code c*N/f} at {@code k/f}, final payment {@code (1 + c/f) N} at maturity.
     * For {@code f = 1} this is exactly {@code C_1 = ... = C_{n-1} = cN, C_n = (1+c)N}.
     */
    public static CashflowSchedule couponBond(String notional, String couponRate, String maturity, int couponFrequency) {
        BigDecimal n = new BigDecimal(notional.trim());
        BigDecimal c = new BigDecimal(couponRate.trim());
        BigDecimal m = new BigDecimal(maturity.trim());
        if (couponFrequency <= 0) {
            throw new IllegalArgumentException("couponFrequency must be positive");
        }
        if (n.signum() <= 0 || c.signum() < 0 || m.signum() <= 0) {
            throw new IllegalArgumentException("notional and maturity must be > 0, coupon rate >= 0");
        }
        BigDecimal f = BigDecimal.valueOf(couponFrequency);
        BigDecimal periods = m.multiply(f);
        if (periods.stripTrailingZeros().scale() > 0) {
            throw new IllegalArgumentException("maturity * couponFrequency must be an integer");
        }
        int count = periods.intValueExact();
        // exact division required: no hidden rounding of contractual amounts
        BigDecimal coupon;
        BigDecimal step;
        try {
            coupon = c.multiply(n).divide(f);
            step = BigDecimal.ONE.divide(f);
        } catch (ArithmeticException e) {
            throw new IllegalArgumentException("coupon amount / period length not exactly representable for frequency "
                    + couponFrequency, e);
        }
        List<BigDecimal> t = new ArrayList<>();
        List<BigDecimal> a = new ArrayList<>();
        for (int k = 1; k <= count; k++) {
            boolean last = k == count;
            BigDecimal amount = last ? coupon.add(n) : coupon;
            if (amount.signum() == 0) {
                continue; // zero coupon rate: only the redemption is a (positive) cashflow
            }
            t.add(last ? m : step.multiply(BigDecimal.valueOf(k)));
            a.add(amount);
        }
        return new CashflowSchedule(t, a);
    }

    /** Single cashflow (zero-coupon bond), e.g. for European reference cases. */
    public static CashflowSchedule single(String time, String amount) {
        return new CashflowSchedule(List.of(new BigDecimal(time.trim())), List.of(new BigDecimal(amount.trim())));
    }

    /** All cashflows scaled by an exact decimal factor (for linear scaling tests and portfolios). */
    public CashflowSchedule scaled(BigDecimal factor) {
        return new CashflowSchedule(Arrays.asList(times), Arrays.stream(amounts).map(x -> x.multiply(factor)).toList());
    }

    public int size() { return times.length; }
    public BigDecimal time(int i) { return times[i]; }
    public BigDecimal amount(int i) { return amounts[i]; }
    public double timeDouble(int i) { return timesDouble[i]; }
    public double amountDouble(int i) { return amountsDouble[i]; }
    public BigDecimal lastTime() { return times[times.length - 1]; }

    /**
     * Index of the first cashflow strictly after {@code time}. Cashflows paid exactly at {@code time}
     * are <em>not</em> part of the remaining stream (binding convention).
     */
    public int firstIndexAfter(BigDecimal time) {
        for (int i = 0; i < times.length; i++) {
            if (times[i].compareTo(time) > 0) {
                return i;
            }
        }
        return times.length;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("CashflowSchedule[");
        for (int i = 0; i < times.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(times[i].toPlainString()).append(':').append(amounts[i].toPlainString());
        }
        return sb.append(']').toString();
    }
}

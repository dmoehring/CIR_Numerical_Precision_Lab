package lab.cir.numeric;

import ch.obermuhlner.math.big.BigDecimalMath;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Elementary functions for the BigDecimal variant, backed by {@code ch.obermuhlner:big-math}.
 *
 * <p>Each function is evaluated with {@code precision + guardDigits} digits and the result is rounded to
 * the working {@link MathContext}. There is deliberately no conversion to {@code double} anywhere in
 * this class.
 */
public final class DecimalMath {

    private final PrecisionSettings settings;
    private final MathContext mc;
    private final MathContext internal;

    public DecimalMath(PrecisionSettings settings) {
        this.settings = settings;
        this.mc = settings.mathContext();
        this.internal = settings.internalMathContext();
    }

    public PrecisionSettings settings() { return settings; }

    public MathContext mc() { return mc; }

    public BigDecimal exp(BigDecimal x) {
        if (x.signum() == 0) {
            return BigDecimal.ONE;
        }
        return BigDecimalMath.exp(x, internal).round(mc);
    }

    /**
     * {@code exp(x) - 1} without catastrophic cancellation for small |x|: the internal precision is
     * increased by the number of leading decimal digits that cancel.
     */
    public BigDecimal expm1(BigDecimal x) {
        if (x.signum() == 0) {
            return BigDecimal.ZERO;
        }
        int magnitude = x.precision() - x.scale() - 1; // floor(log10|x|)
        int extra = Math.max(0, -magnitude);
        MathContext wide = new MathContext(internal.getPrecision() + extra, internal.getRoundingMode());
        return BigDecimalMath.exp(x, wide).subtract(BigDecimal.ONE, wide).round(mc);
    }

    public BigDecimal log(BigDecimal x) {
        if (x.signum() <= 0) {
            throw new ArithmeticException("log of non-positive value " + x);
        }
        return BigDecimalMath.log(x, internal).round(mc);
    }

    public BigDecimal sqrt(BigDecimal x) {
        int s = x.signum();
        if (s < 0) {
            throw new ArithmeticException("sqrt of negative value " + x);
        }
        if (s == 0) {
            return BigDecimal.ZERO;
        }
        return BigDecimalMath.sqrt(x, internal).round(mc);
    }

    /** {@code max(x, 0)} (no rounding involved). */
    public static BigDecimal positivePart(BigDecimal x) {
        return x.signum() < 0 ? BigDecimal.ZERO : x;
    }
}

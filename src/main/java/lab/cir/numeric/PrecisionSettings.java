package lab.cir.numeric;

import java.math.MathContext;
import java.math.RoundingMode;

/**
 * Decimal precision strategy.
 *
 * <ul>
 *   <li><b>working precision</b> ({@link #mathContext()}): applied to every arithmetic operation of the
 *       BigDecimal valuation (add, subtract, multiply, divide) and to the results of elementary functions;</li>
 *   <li><b>guard digits</b>: elementary functions (exp, log, sqrt) are evaluated internally with
 *       {@code precision + guardDigits} digits and only then rounded to the working precision, so their
 *       results are (almost always) correctly rounded at working precision;</li>
 *   <li><b>output rounding</b> is not done here: reports print the full decimal value.</li>
 * </ul>
 */
public record PrecisionSettings(int precision, RoundingMode roundingMode, int guardDigits) {

    public static final PrecisionSettings DECIMAL128 = new PrecisionSettings(34, RoundingMode.HALF_EVEN, 10);

    public PrecisionSettings {
        if (precision < 1) {
            throw new IllegalArgumentException("precision must be >= 1 (unlimited precision is not supported)");
        }
        if (guardDigits < 0) {
            throw new IllegalArgumentException("guardDigits must be >= 0");
        }
    }

    public MathContext mathContext() {
        return new MathContext(precision, roundingMode);
    }

    public MathContext internalMathContext() {
        return new MathContext(precision + guardDigits, roundingMode);
    }

    public PrecisionSettings withPrecision(int newPrecision) {
        return new PrecisionSettings(newPrecision, roundingMode, guardDigits);
    }

    @Override
    public String toString() {
        return precision + " digits, " + roundingMode + ", guard " + guardDigits;
    }
}

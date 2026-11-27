package lab.cir.analytic;

/**
 * Algebraically equivalent evaluation forms of {@code ln A(tau)}; they differ only in floating-point behaviour.
 * With {@code g = (theta + h) tau / 2}, {@code em1 = e^{h tau} - 1}, {@code denom = 2h + (theta + h) em1}:
 * <ul>
 *   <li>{@link #LOG_DIFFERENCE} (textbook, default): {@code ln A = c (ln 2h + g - ln denom)}. The bracket is a small
 *       difference of O(1) terms; in double this cancellation costs up to ~1e-15 relative accuracy of P for
 *       maturities of a few years (measured, see docs).</li>
 *   <li>{@link #LOG1P_RATIO}: {@code ln A = c log1p((2h expm1(g) - (theta + h) em1) / denom)}, i.e. the logarithm of
 *       the ratio {@code 2h e^g / denom} close to 1 evaluated through its distance to 1. Better for short and medium
 *       maturities; for very long maturities the numerator itself cancels.</li>
 * </ul>
 */
public enum ZeroBondForm {
    LOG_DIFFERENCE,
    LOG1P_RATIO
}

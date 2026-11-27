package lab.cir.analytic;

import lab.cir.model.CirParameters;

/**
 * Analytic CIR zero-coupon bond price in {@code double}:
 *
 * <pre>
 *   P(t,T) = A(tau) exp(-B(tau) r_t),   tau = T - t,   h = sqrt(theta^2 + 2 sigma^2)
 *   B(tau) = 2 (e^{h tau} - 1) / (2h + (theta + h)(e^{h tau} - 1))
 *   A(tau) = [ 2h e^{(theta + h) tau / 2} / (2h + (theta + h)(e^{h tau} - 1)) ]^{2 theta mu / sigma^2}
 * </pre>
 *
 * Evaluated as {@code exp(logA - B r)} with {@code logA = (2 theta mu / sigma^2) (ln 2h + (theta+h) tau/2 - ln denom)}
 * and {@code expm1} for {@code e^{h tau} - 1}, which avoids overflow of the power for large exponents
 * {@code 2 theta mu / sigma^2} (small volatility) and cancellation for short maturities.
 */
public final class CirZeroBondDouble {

    /** {@code ln A(tau)} and {@code B(tau)} for a fixed time to maturity. */
    public record Coefficients(double logA, double b) {
        public double price(double rate) {
            return Math.exp(logA - b * rate);
        }
    }

    private final double theta;
    private final double sigma;
    private final double h;
    private final double exponent;
    private final ZeroBondForm form;

    public CirZeroBondDouble(CirParameters p) {
        this(p, ZeroBondForm.LOG_DIFFERENCE);
    }

    public CirZeroBondDouble(CirParameters p, ZeroBondForm form) {
        this.form = form;
        this.theta = p.meanReversionSpeedDouble();
        double mu = p.longTermMeanDouble();
        this.sigma = p.volatilityDouble();
        this.h = Math.sqrt(theta * theta + 2.0 * sigma * sigma);
        this.exponent = 2.0 * theta * mu / (sigma * sigma);
    }

    public Coefficients coefficients(double tau) {
        if (tau < 0) {
            throw new IllegalArgumentException("tau must be >= 0");
        }
        if (tau == 0) {
            return new Coefficients(0.0, 0.0);
        }
        double em1 = Math.expm1(h * tau);
        double denom = 2.0 * h + (theta + h) * em1;
        double b = 2.0 * em1 / denom;
        double logA;
        if (form == ZeroBondForm.LOG_DIFFERENCE) {
            logA = exponent * (Math.log(2.0 * h) + (theta + h) * tau / 2.0 - Math.log(denom));
        } else {
            double g = (theta + h) * tau / 2.0;
            logA = exponent * Math.log1p((2.0 * h * Math.expm1(g) - (theta + h) * em1) / denom);
        }
        return new Coefficients(logA, b);
    }

    public double price(double tau, double rate) {
        return coefficients(tau).price(rate);
    }

    public double h() { return h; }

    public ZeroBondForm form() { return form; }
}

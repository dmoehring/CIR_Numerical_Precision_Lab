package lab.cir.analytic;

import org.apache.commons.math3.special.Gamma;

/**
 * Noncentral chi-square CDF in <b>double</b> (reference implementation, not a high-precision one):
 *
 * <pre>
 *   F(x; k, lambda) = sum_{i >= 0} e^{-lambda/2} (lambda/2)^i / i!  *  P(k/2 + i, x/2)
 * </pre>
 *
 * with the regularized lower incomplete gamma function {@code P} from Apache Commons Math. The Poisson
 * mixture is summed outwards from its mode {@code floor(lambda/2)} (weights computed in log space at the mode and
 * by recurrence elsewhere), which avoids underflow of {@code e^{-lambda/2}} for large noncentrality.
 * Summation stops when the Poisson weights fall below {@code 1e-17}; the sum is normalised by the Poisson mass used.
 * Validated against SciPy {@code ncx2.cdf} (see {@code validation/reference_values.py}).
 */
public final class NoncentralChiSquared {

    private static final double GAMMA_EPS = 1e-16;
    private static final int GAMMA_MAX_ITER = 100_000;
    private static final double TAIL = 1e-17;

    private NoncentralChiSquared() {
    }

    public static double cdf(double x, double degreesOfFreedom, double noncentrality) {
        if (x <= 0.0) {
            return 0.0;
        }
        if (degreesOfFreedom <= 0.0 || noncentrality < 0.0) {
            throw new IllegalArgumentException("invalid parameters");
        }
        double halfK = degreesOfFreedom / 2.0;
        double halfX = x / 2.0;
        if (noncentrality == 0.0) {
            return Gamma.regularizedGammaP(halfK, halfX, GAMMA_EPS, GAMMA_MAX_ITER);
        }
        double mean = noncentrality / 2.0;
        int mode = (int) Math.floor(mean);
        double logW = -mean + mode * Math.log(mean) - Gamma.logGamma(mode + 1.0);
        double wMode = Math.exp(logW);

        double sum = 0.0;
        double massUsed = 0.0;
        // forward from the mode
        double w = wMode;
        for (int i = mode; ; i++) {
            sum += w * Gamma.regularizedGammaP(halfK + i, halfX, GAMMA_EPS, GAMMA_MAX_ITER);
            massUsed += w;
            w *= mean / (i + 1.0);
            if (w < TAIL && i > mode + 10) {
                break;
            }
            if (i - mode > 10_000_000) {
                throw new ArithmeticException("noncentral chi-square series did not converge");
            }
        }
        // backward from the mode
        w = wMode;
        for (int i = mode - 1; i >= 0; i--) {
            w *= (i + 1.0) / mean;
            sum += w * Gamma.regularizedGammaP(halfK + i, halfX, GAMMA_EPS, GAMMA_MAX_ITER);
            massUsed += w;
            if (w < TAIL && i < mode - 10) {
                break;
            }
        }
        if (Math.abs(massUsed - 1.0) > 1e-9) {
            throw new ArithmeticException("Poisson weights do not sum to 1 (" + massUsed + ")");
        }
        // normalise by the Poisson mass actually used: removes the common relative error of the anchor weight
        // exp(logW) (log-gamma of a large mode has an absolute error of ~ulp(lgamma), i.e. ~1e-11 relative for lambda ~ 1e4)
        return Math.min(1.0, sum / massUsed);
    }
}

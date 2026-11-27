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

    // ---------------------------------------------------------------------------------------------------------
    // density

    /**
     * Noncentral chi-square density {@code f(x; k, lambda) = sum_i Pois(i; lambda/2) g_{k+2i}(x)} with the central
     * densities {@code g_n(x) = Gamma(x/2; n/2) / 2}. Every term is a product of two Poisson-type densities
     * evaluated with Loader's saddle-point form ({@code exp(-stirlerr(n) - bd0(n, m)) / sqrt(2 pi n)}), which is
     * accurate to a few ulp even for large arguments (no cancellation between {@code lambda/2}, {@code i ln(lambda/2)}
     * and {@code ln i!}). The series is summed outwards from its largest term until terms drop below
     * {@code 1e-17} of the sum.
     */
    public static double density(double x, double degreesOfFreedom, double noncentrality) {
        if (degreesOfFreedom <= 0.0 || noncentrality < 0.0) {
            throw new IllegalArgumentException("invalid parameters");
        }
        if (x < 0.0) {
            return 0.0;
        }
        double halfK = degreesOfFreedom / 2.0;
        double mean = noncentrality / 2.0;
        double halfX = x / 2.0;
        if (x == 0.0) {
            if (halfK > 1.0) return 0.0;
            if (halfK == 1.0) return 0.5 * Math.exp(-mean);
            return Double.POSITIVE_INFINITY;
        }
        if (mean == 0.0) {
            return 0.5 * gammaDensity(halfK, halfX);
        }
        // index of the largest term: ratio t_{i+1}/t_i = mean/(i+1) * halfX/(halfK+i) = 1
        double b = halfK + 1.0;
        double disc = b * b - 4.0 * (halfK - mean * halfX);
        int peak = disc > 0 ? (int) Math.max(0.0, Math.round((-b + Math.sqrt(disc)) / 2.0)) : 0;
        double tPeak = poisson(peak, mean) * 0.5 * gammaDensity(halfK + peak, halfX);
        if (tPeak == 0.0) {
            return 0.0;
        }
        double sum = 1.0;
        double t = 1.0;
        for (int i = peak; ; i++) {
            t *= mean / (i + 1.0) * halfX / (halfK + i);
            sum += t;
            if (t < TAIL * sum) break;
        }
        t = 1.0;
        for (int i = peak; i > 0; i--) {
            t *= i / mean * (halfK + i - 1.0) / halfX;
            sum += t;
            if (t < TAIL * sum) break;
        }
        return tPeak * sum;
    }

    /** Gamma(shape, 1) density at {@code y > 0} (Loader). */
    static double gammaDensity(double shape, double y) {
        if (shape >= 1.0) {
            return poisson(shape - 1.0, y);
        }
        return poisson(shape, y) * shape / y;
    }

    /** {@code m^n e^{-m} / Gamma(n+1)} for real {@code n >= 0}, Loader's saddle-point form. */
    static double poisson(double n, double m) {
        if (m == 0.0) return n == 0.0 ? 1.0 : 0.0;
        if (n == 0.0) return Math.exp(-m);
        return Math.exp(-stirlerr(n) - bd0(n, m)) / Math.sqrt(2.0 * Math.PI * n);
    }

    /** {@code n ln(n/m) + m - n}, computed without cancellation for {@code n ~ m}. */
    static double bd0(double n, double m) {
        if (Math.abs(n - m) < 0.1 * (n + m)) {
            double v = (n - m) / (n + m);
            double s = (n - m) * v;
            double ej = 2.0 * n * v;
            v *= v;
            for (int j = 1; j < 1000; j++) {
                ej *= v;
                double s1 = s + ej / (2 * j + 1);
                if (s1 == s) return s1;
                s = s1;
            }
            return s;
        }
        return n * Math.log(n / m) + m - n;
    }

    /** Stirling error {@code ln Gamma(n+1) - (n + 1/2) ln n + n - ln sqrt(2 pi)}. */
    static double stirlerr(double n) {
        final double s0 = 1.0 / 12, s1 = 1.0 / 360, s2 = 1.0 / 1260, s3 = 1.0 / 1680, s4 = 1.0 / 1188;
        if (n <= 15.0) {
            double nn = n + n;
            if (nn == Math.rint(nn)) return STIRLERR_HALVES[(int) nn];
            return Gamma.logGamma(n + 1.0) - (n + 0.5) * Math.log(n) + n - 0.5 * Math.log(2.0 * Math.PI);
        }
        double nn = n * n;
        if (n > 500) return (s0 - s1 / nn) / n;
        if (n > 80) return (s0 - (s1 - s2 / nn) / nn) / n;
        if (n > 35) return (s0 - (s1 - (s2 - s3 / nn) / nn) / nn) / n;
        return (s0 - (s1 - (s2 - (s3 - s4 / nn) / nn) / nn) / nn) / n;
    }

    /** stirlerr(k/2), k = 0..30, computed with mpmath (40 digits). */
    private static final double[] STIRLERR_HALVES = {
            0,
            0.1534264097200273452914, 0.08106146679532725821967, 0.05481412105191765389614,
            0.04134069595540929409382, 0.03316287351993628748511, 0.02767792568499833914879,
            0.02374616365629749597133, 0.02079067210376509311152, 0.01848845053267318523078,
            0.01664469118982119216319, 0.01513497322191737887351, 0.01387612882307074799875,
            0.01281046524292022692425, 0.01189670994589177009506, 0.01110455975820691732663,
            0.01041126526197209649748, 0.00979941612615880329839, 0.009255462182712732917729,
            0.008768700134139385462955, 0.008330563433362871256469, 0.00793411456431402054725,
            0.007573675487951840794972, 0.007244554301320383179546, 0.006942840107209529865664,
            0.006665247032707682442356, 0.00640899418800420706844, 0.006171712263039457647535,
            0.005951370112758847735624, 0.005746216513010115682026, 0.005554733551962801371039,
    };
}

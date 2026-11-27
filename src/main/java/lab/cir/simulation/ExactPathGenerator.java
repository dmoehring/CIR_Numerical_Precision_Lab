package lab.cir.simulation;

import java.util.SplittableRandom;

import lab.cir.analytic.CirZeroBondDouble;
import lab.cir.analytic.RemainingBondDouble;
import lab.cir.analytic.ZeroBondForm;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;
import lab.cir.random.InverseNormal;
import lab.cir.random.NormalSource;

/**
 * Exact simulation of the CIR short rate at the exercise dates (no time discretisation), under the forward
 * measure {@code Q^S} of the <b>last exercise date</b> {@code S = e_m}, whose numeraire is the zero bond
 * {@code P(t, S)}.
 *
 * <p>Under {@code Q^S} the transition {@code r_t -> r_T} ({@code t < T <= S}, {@code tau = T - t}) is a scaled
 * noncentral chi-square (Brigo &amp; Mercurio 2006, eq. (3.28)):
 * <pre>
 *   2 (rho + psi + B(S - T)) r_T  ~  chi'^2( 4 theta mu / sigma^2,  2 rho^2 r_t e^{h tau} / (rho + psi + B(S - T)) )
 *   rho = 2h / (sigma^2 (e^{h tau} - 1)),   psi = (theta + h) / sigma^2
 * </pre>
 * sampled as a Poisson mixture: {@code N ~ Poisson(lambda/2)}, {@code chi'^2 = 2 Gamma(k/2 + N)}
 * (Marsaglia-Tsang gamma sampler with AS241 normals, inversion for the Poisson variable).
 *
 * <p>Change of numeraire: a cashflow {@code X} at {@code e_j} is worth {@code P(e_i, S) E^S[X / P(e_j, S)]} at
 * {@code e_i}. The generator therefore stores {@code I_j = ln(P(e_j, S; r_j) / P(0, S; r_0))} as
 * {@link PathDataDouble#cumulativeIntegral}: the induction's discount factor {@code exp(-(I_tau - I_j))} becomes
 * {@code P(e_j, S) / P(tau, S)} and {@code exp(-I_tau)} becomes {@code P(0, S) / P(tau, S)}, so the existing
 * backward induction and policy valuation work unchanged and estimate the same continuation values
 * {@code C_j(r)} as under the risk-neutral measure. Only the distribution of the simulated states differs (it is
 * the {@code Q^S} law, slightly shifted towards lower rates).
 *
 * <p>Random numbers: one {@link SplittableRandom} per path seeded like {@link NormalSource}, so results do not
 * depend on the thread count.
 */
public final class ExactPathGenerator {

    private final double r0;
    private final double dof;
    private final BermudanCallOption option;
    private final RemainingBondDouble bond;
    private final double[] scale;
    private final double[] nonCentralityPerRate;
    private final CirZeroBondDouble.Coefficients[] toLast;
    private final double logP0;

    public ExactPathGenerator(CirParameters p, BermudanCallOption option) {
        this.option = option;
        double theta = p.meanReversionSpeedDouble();
        double mu = p.longTermMeanDouble();
        double sigma = p.volatilityDouble();
        this.r0 = p.initialRateDouble();
        double s2 = sigma * sigma;
        this.dof = 4.0 * theta * mu / s2;
        CirZeroBondDouble zb = new CirZeroBondDouble(p, ZeroBondForm.LOG1P_RATIO);
        double h = zb.h();
        this.bond = new RemainingBondDouble(option, zb);
        int m = option.exerciseDateCount();
        double last = option.exerciseDateDouble(m - 1);
        this.scale = new double[m];
        this.nonCentralityPerRate = new double[m];
        this.toLast = new CirZeroBondDouble.Coefficients[m];
        double psi = (theta + h) / s2;
        for (int j = 0; j < m; j++) {
            double from = j == 0 ? 0.0 : option.exerciseDateDouble(j - 1);
            double to = option.exerciseDateDouble(j);
            double tau = to - from;
            double rho = 2.0 * h / (s2 * Math.expm1(h * tau));
            double bRest = zb.coefficients(last - to).b();
            scale[j] = 2.0 * (rho + psi + bRest);
            nonCentralityPerRate[j] = 2.0 * rho * rho * Math.exp(h * tau) / (rho + psi + bRest);
            toLast[j] = zb.coefficients(last - to);
        }
        CirZeroBondDouble.Coefficients c0 = zb.coefficients(last);
        this.logP0 = c0.logA() - c0.b() * r0;
    }

    /** Rates, numeraire-based discount exponents, bond and exercise values at the exercise dates. */
    public PathDataDouble simulate(long seed, int paths, int threads) {
        int m = option.exerciseDateCount();
        PathDataDouble d = new PathDataDouble(m, paths);
        ParallelPaths.forEach(paths, threads, p -> {
            SplittableRandom rng = new SplittableRandom(NormalSource.pathSeed(seed, p));
            double r = r0;
            for (int j = 0; j < m; j++) {
                double lambda = nonCentralityPerRate[j] * r;
                r = noncentralChiSquare(rng, dof, lambda) / scale[j];
                d.rate[j][p] = r;
                double logP = toLast[j].logA() - toLast[j].b() * r;
                d.cumulativeIntegral[j][p] = logP - logP0;
                double b = bond.bondValue(j, r);
                d.bondValue[j][p] = b;
                d.exerciseValue[j][p] = bond.exerciseValue(b);
            }
        });
        return d;
    }

    /** {@code P(0, S)}: converts expectations of {@code X / P(e_j, S)} into prices. */
    public double numeraireAtZero() {
        return Math.exp(logP0);
    }

    // ------------------------------------------------------------------------------------------ samplers

    static double noncentralChiSquare(SplittableRandom rng, double k, double lambda) {
        int n = lambda > 0.0 ? poisson(rng, 0.5 * lambda) : 0;
        return 2.0 * gamma(rng, 0.5 * k + n);
    }

    /** Poisson by sequential inversion (means here are below ~100). */
    static int poisson(SplittableRandom rng, double mean) {
        if (mean > 600.0) {
            throw new IllegalStateException("Poisson mean " + mean + " outside the supported range");
        }
        double u = NormalSource.uniform(rng.nextLong());
        double p = Math.exp(-mean);
        double cdf = p;
        int n = 0;
        while (u > cdf) {
            n++;
            p *= mean / n;
            cdf += p;
            if (p < 1e-300 && n > mean) break;
        }
        return n;
    }

    /** Gamma(shape, 1) by Marsaglia and Tsang (2000); shape < 1 via {@code Gamma(a+1) U^{1/a}}. */
    static double gamma(SplittableRandom rng, double shape) {
        if (shape < 1.0) {
            double u = NormalSource.uniform(rng.nextLong());
            return gamma(rng, shape + 1.0) * Math.pow(u, 1.0 / shape);
        }
        double d = shape - 1.0 / 3.0;
        double c = 1.0 / Math.sqrt(9.0 * d);
        while (true) {
            double x = InverseNormal.quantile(NormalSource.uniform(rng.nextLong()));
            double v = 1.0 + c * x;
            if (v <= 0.0) continue;
            v = v * v * v;
            double u = NormalSource.uniform(rng.nextLong());
            if (Math.log(u) < 0.5 * x * x + d - d * v + d * Math.log(v)) {
                return d * v;
            }
        }
    }
}

package lab.cir.analytic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;

/**
 * Semi-analytic double-precision reference for the Bermudan bond call, independent of the finite-difference solver:
 * backward induction over the exercise dates only, with the continuation value computed by numerical integration
 * against the <b>exact</b> CIR transition law.
 *
 * <p>Between two exercise dates {@code e_j < e_{j+1}} (length {@code tau}) the continuation value is written under the
 * {@code e_{j+1}}-forward measure, whose numeraire is the zero bond maturing at {@code e_{j+1}}:
 *
 * <pre>
 *   C_j(r) = E^Q[ exp(-int r ds) V_{j+1}(r_{e_{j+1}}) | r_{e_j} = r ] = P(tau; r) * E^{T}[ V_{j+1}(r_{e_{j+1}}) | r ]
 * </pre>
 *
 * so no path integral is needed. Under this measure {@code 2 (rho + psi) r_{e_{j+1}}} is noncentral chi-square with
 * {@code 4 theta mu / sigma^2} degrees of freedom and noncentrality {@code 2 rho^2 r e^{h tau} / (rho + psi)},
 * {@code rho = 2h / (sigma^2 (e^{h tau} - 1))}, {@code psi = (theta + h) / sigma^2} (Cox, Ingersoll, Ross 1985;
 * Brigo &amp; Mercurio 2006, Sect. 3.2.3; the same law yields the closed-form zero-bond option price).
 *
 * <ul>
 *   <li>{@code C_j} is stored as a piecewise Chebyshev interpolant (equal panels on {@code [0, rateMax]},
 *       {@code degree} first-kind Chebyshev nodes per panel, barycentric evaluation). It is analytic in {@code r},
 *       so the interpolation converges geometrically; local panels keep the (irrelevant) truncation effects at
 *       {@code rateMax} local. Beyond {@code rateMax} the last value is held constant (weight below 1e-15).</li>
 *   <li>{@code V_{j+1} = max(C_{j+1}, E_{j+1})} is integrated piecewise: the kinks ({@code C = E}, located by bisection
 *       on the interpolant and the exact exercise value, and {@code B = K}) and all panel boundaries are break
 *       points, so every piece has a smooth integrand. Each piece uses adaptive Gauss-Kronrod (7/15).</li>
 *   <li>The density is the Poisson mixture of {@link NoncentralChiSquared#density}; zero bonds use the
 *       {@link ZeroBondForm#LOG1P_RATIO} form.</li>
 * </ul>
 * There is no time-stepping and no Euler bias; the remaining errors (interpolation, quadrature, truncation, double
 * rounding) are controlled by {@link Settings} and estimated by varying them.
 *
 * <p><b>Policy valuation.</b> {@link #valuePolicy} values an arbitrary exercise rule (e.g. one trained by quantile
 * clusters or least squares) exactly in the continuous-time model: {@code V^pi_j = E_j} where the rule exercises and
 * {@code C^pi_j = P E^T[V^pi_{j+1}]} elsewhere. The switch points of the rule are located by scanning its decisions on
 * a fine grid and bisection, and become break points of the integration. The optimal price is the special case of
 * the rule "exercise iff {@code E > C}". Instances are thread-safe.
 */
public final class CirQuadratureBermudan {

    /**
     * @param rateMax   upper end of the interpolation range for the continuation values
     * @param panels    number of equal Chebyshev panels on {@code [0, rateMax]}
     * @param degree    Chebyshev nodes per panel
     * @param tolerance absolute local error target of the adaptive quadrature (per unit of the integration range)
     */
    public record Settings(double rateMax, int panels, int degree, double tolerance) {
        public static Settings standard() {
            return new Settings(1.0, 64, 16, 1e-13);
        }
    }

    /** Exercise rule at the exercise dates (date index, rate, remaining bond value, exercise value). */
    @FunctionalInterface
    public interface ExercisePolicy {
        boolean exercise(int date, double rate, double bondValue, double exerciseValue);
    }

    public final class Result {
        private final double price;
        private final double[][] boundaries;
        private final Chebyshev[] continuation;
        private final long densityEvaluations;
        private final Settings settings;

        private Result(double price, double[][] boundaries, Chebyshev[] continuation, long densityEvaluations,
                       Settings settings) {
            this.price = price;
            this.boundaries = boundaries;
            this.continuation = continuation;
            this.densityEvaluations = densityEvaluations;
            this.settings = settings;
        }

        public double price() { return price; }

        /** Per exercise date the switch points of the rule (optimal rule: {@code C_j = E_j}; last date: {@code B = K}). */
        public double[][] exerciseBoundaries() { return boundaries; }

        /** Continuation value {@code C_j(r)} (of the valued rule) at exercise date {@code j < m-1} (interpolant). */
        public double continuationValue(int j, double r) { return continuation[j].value(r); }

        public long densityEvaluations() { return densityEvaluations; }

        public Settings settings() { return settings; }
    }

    private final double theta;
    private final double mu;
    private final double sigma;
    private final double r0;
    private final double h;
    private final CirZeroBondDouble zeroBond;

    public CirQuadratureBermudan(CirParameters p) {
        this.theta = p.meanReversionSpeedDouble();
        this.mu = p.longTermMeanDouble();
        this.sigma = p.volatilityDouble();
        this.r0 = p.initialRateDouble();
        this.zeroBond = new CirZeroBondDouble(p, ZeroBondForm.LOG1P_RATIO);
        this.h = zeroBond.h();
    }

    public Result price(BermudanCallOption option) {
        return price(option, Settings.standard());
    }

    public Result price(BermudanCallOption option, Settings s) {
        return valuePolicy(option, null, s);
    }

    /**
     * Value of an exercise rule in the continuous-time model ({@code policy == null}: optimal rule, i.e. the price).
     * The rule is evaluated at rates in {@code [0, rateMax]}; beyond, its decision at {@code rateMax} is kept.
     */
    public Result valuePolicy(BermudanCallOption option, ExercisePolicy policy, Settings s) {
        long[] evaluations = new long[1];
        RemainingBondDouble bond = new RemainingBondDouble(option, zeroBond);
        int m = option.exerciseDateCount();
        double[][] boundaries = new double[m][];
        Chebyshev[] cont = new Chebyshev[m];

        ValueFunction next = new ValueFunction(m - 1, null, bond, s, policy);
        boundaries[m - 1] = next.switchPoints();
        for (int j = m - 2; j >= 0; j--) {
            Transition tr = new Transition(option.exerciseDateDouble(j + 1) - option.exerciseDateDouble(j), evaluations);
            Chebyshev cj = new Chebyshev(s.rateMax(), s.panels(), s.degree());
            final ValueFunction vf = next;
            cj.fill(r -> tr.discount(r) * integrate(vf, tr, r, s));
            cont[j] = cj;
            next = new ValueFunction(j, cj, bond, s, policy);
            boundaries[j] = next.switchPoints();
        }
        Transition first = new Transition(option.exerciseDateDouble(0), evaluations);
        double price = first.discount(r0) * integrate(next, first, r0, s);
        return new Result(price, boundaries, cont, evaluations[0], s);
    }

    // ------------------------------------------------------------------------------------------ transition law

    /** Law of {@code r_{t+tau}} given {@code r_t} under the {@code (t+tau)}-forward measure. */
    final class Transition {
        final double scale;
        final double dof;
        final double nonCentralityPerRate;
        final CirZeroBondDouble.Coefficients bondCoefficients;
        private final long[] evaluations;

        Transition(double tau, long[] evaluations) {
            this.evaluations = evaluations;
            double s2 = sigma * sigma;
            double rho = 2.0 * h / (s2 * Math.expm1(h * tau));
            double psi = (theta + h) / s2;
            this.scale = 2.0 * (rho + psi);
            this.dof = 4.0 * theta * mu / s2;
            this.nonCentralityPerRate = 2.0 * rho * rho * Math.exp(h * tau) / (rho + psi);
            this.bondCoefficients = zeroBond.coefficients(tau);
        }

        double discount(double r) {
            return bondCoefficients.price(r);
        }

        double density(double x, double r) {
            evaluations[0]++;
            return scale * NoncentralChiSquared.density(scale * x, dof, nonCentralityPerRate * r);
        }

        /** Integration range in the rate variable: all but a negligible part of the probability mass. */
        double[] support(double r) {
            double lambda = nonCentralityPerRate * r;
            double meanY = dof + lambda;
            double sdY = Math.sqrt(2.0 * (dof + 2.0 * lambda));
            double lo = Math.max(0.0, meanY - 25.0 * sdY);
            double hi = meanY + 40.0 * sdY + 60.0;
            return new double[]{lo / scale, hi / scale};
        }
    }

    // ------------------------------------------------------------------------------------------ value function

    /**
     * {@code V_j = E_j} where the rule exercises, else {@code C_j} (0 at the last date); with the optimal rule
     * {@code V_j = max(C_j, E_j)}. Break points: switch points of the rule and the kink {@code B = K} of {@code E_j}.
     */
    private static final class ValueFunction {
        private static final int MIN_SCAN_POINTS = 20_000;

        private final int date;
        private final Chebyshev continuation;
        private final RemainingBondDouble bond;
        private final ExercisePolicy policy;
        private final double[] switchPoints;
        private final double[] kinks;
        private final double rateMax;

        ValueFunction(int date, Chebyshev continuation, RemainingBondDouble bond, Settings s, ExercisePolicy rule) {
            this.date = date;
            this.continuation = continuation;
            this.bond = bond;
            this.rateMax = s.rateMax();
            this.policy = rule != null ? rule
                    : (j, r, b, e) -> e > 0.0 && e > (continuation == null ? 0.0 : continuation.value(r));
            // switch points: decision changes on a fine scan of [0, rateMax], refined by bisection
            List<Double> sw = new ArrayList<>();
            int scan = Math.max(MIN_SCAN_POINTS, 4 * s.degree() * s.panels());
            double prevX = 0.0;
            boolean prevD = decide(0.0);
            for (int k = 1; k <= scan; k++) {
                double x = s.rateMax() * k / scan;
                boolean d = decide(x);
                if (d != prevD) {
                    final boolean left = prevD;
                    sw.add(bisect(t -> decide(t) == left ? 1.0 : -1.0, prevX, x));
                }
                prevX = x;
                prevD = d;
            }
            this.switchPoints = sw.stream().mapToDouble(Double::doubleValue).toArray();
            List<Double> all = new ArrayList<>(sw);
            double strikeKink = strikeKink(s.rateMax());
            if (!Double.isNaN(strikeKink)) all.add(strikeKink);
            this.kinks = all.stream().mapToDouble(Double::doubleValue).sorted().toArray();
        }

        private boolean decide(double x) {
            double b = bond.bondValue(date, x);
            return policy.exercise(date, x, b, bond.exerciseValue(b));
        }

        double value(double x) {
            double xr = Math.min(x, rateMax);
            double b = bond.bondValue(date, x);
            double e = bond.exerciseValue(b);
            double bDecision = xr == x ? b : bond.bondValue(date, xr);
            if (policy.exercise(date, xr, bDecision, bond.exerciseValue(bDecision))) {
                return e;
            }
            return continuation == null ? 0.0 : continuation.value(xr);
        }

        double[] kinks() { return kinks; }

        double[] switchPoints() { return switchPoints; }

        double[] panelBreaks() { return continuation == null ? new double[0] : continuation.breaks(); }

        private double strikeKink(double rMax) {
            double k = bond.strike();
            if (bond.bondValue(date, 0.0) <= k) return Double.NaN;
            double hi = rMax;
            while (bond.bondValue(date, hi) > k) hi *= 2.0;
            return bisect(x -> bond.bondValue(date, x) - k, 0.0, hi);
        }
    }

    // ------------------------------------------------------------------------------------------ integration

    private double integrate(ValueFunction v, Transition tr, double r, Settings s) {
        double[] range = tr.support(r);
        double lo = range[0];
        double hi = range[1];
        TreeSet<Double> breaks = new TreeSet<>();
        breaks.add(lo);
        breaks.add(hi);
        for (double x : v.kinks()) if (x > lo && x < hi) breaks.add(x);
        for (double x : v.panelBreaks()) if (x > lo && x < hi) breaks.add(x);
        double span = hi - lo;
        double total = 0.0;
        Double prev = null;
        for (double x : breaks) {
            if (prev != null) {
                final double a = prev;
                total += adaptive(t -> v.value(t) * tr.density(t, r), a, x, s.tolerance() * (x - a), span, 0);
            }
            prev = x;
        }
        return total;
    }

    @FunctionalInterface
    interface Fn {
        double at(double x);
    }

    private static final double[] XGK = {
            0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
            0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
            0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
            0.207784955007898467600689403773245, 0.0};
    private static final double[] WGK = {
            0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
            0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
            0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
            0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
    private static final double[] WG = {
            0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
            0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

    /** Adaptive Gauss-Kronrod 7/15 with bisection until {@code |K15 - G7| <= tol}. */
    static double adaptive(Fn f, double a, double b, double tol, double span, int depth) {
        double c = 0.5 * (a + b);
        double hw = 0.5 * (b - a);
        double fc = f.at(c);
        double k = WGK[7] * fc;
        double g = WG[3] * fc;
        for (int i = 0; i < 7; i++) {
            double dx = hw * XGK[i];
            double f1 = f.at(c - dx);
            double f2 = f.at(c + dx);
            k += WGK[i] * (f1 + f2);
            if (i % 2 == 1) g += WG[i / 2] * (f1 + f2);
        }
        k *= hw;
        g *= hw;
        double err = Math.abs(k - g);
        if (err <= Math.max(tol, 1e-15 * Math.abs(k)) || depth >= 40) {
            return k;
        }
        return adaptive(f, a, c, 0.5 * tol, span, depth + 1) + adaptive(f, c, b, 0.5 * tol, span, depth + 1);
    }

    static double bisect(Fn g, double a, double b) {
        double ga = g.at(a);
        for (int it = 0; it < 200; it++) {
            double mid = 0.5 * (a + b);
            if (mid <= a || mid >= b) break;
            double gm = g.at(mid);
            if ((gm > 0) == (ga > 0)) {
                a = mid;
                ga = gm;
            } else {
                b = mid;
            }
        }
        return 0.5 * (a + b);
    }

    // ------------------------------------------------------------------------------------------ interpolation

    /** Piecewise Chebyshev interpolant on equal panels of {@code [0, rMax]} (first-kind nodes, barycentric form). */
    static final class Chebyshev {
        private final double rMax;
        private final int panels;
        private final double width;
        private final double[] t;
        private final double[] w;
        private final double[][] values;

        Chebyshev(double rMax, int panels, int degree) {
            this.rMax = rMax;
            this.panels = panels;
            this.width = rMax / panels;
            this.t = new double[degree];
            this.w = new double[degree];
            for (int k = 0; k < degree; k++) {
                double angle = (2.0 * k + 1.0) * Math.PI / (2.0 * degree);
                t[k] = Math.cos(angle);
                w[k] = ((k & 1) == 0 ? 1.0 : -1.0) * Math.sin(angle);
            }
            this.values = new double[panels][degree];
        }

        void fill(Fn f) {
            for (int p = 0; p < panels; p++) {
                double mid = (p + 0.5) * width;
                for (int k = 0; k < t.length; k++) {
                    values[p][k] = f.at(mid + 0.5 * width * t[k]);
                }
            }
        }

        double value(double x) {
            if (x <= 0.0) x = 0.0;
            if (x >= rMax) x = rMax;
            int p = Math.min(panels - 1, (int) (x / width));
            double u = (x - (p + 0.5) * width) / (0.5 * width);
            double num = 0.0;
            double den = 0.0;
            double[] y = values[p];
            for (int k = 0; k < t.length; k++) {
                double d = u - t[k];
                if (d == 0.0) return y[k];
                double q = w[k] / d;
                num += q * y[k];
                den += q;
            }
            return num / den;
        }

        double[] breaks() {
            double[] b = new double[panels + 1];
            for (int p = 0; p <= panels; p++) b[p] = p * width;
            return b;
        }

        @Override
        public String toString() {
            return "Chebyshev{panels=" + panels + ", degree=" + t.length + ", values=" + Arrays.deepToString(values) + "}";
        }
    }
}

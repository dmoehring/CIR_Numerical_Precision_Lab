package lab.cir.analytic;

import java.util.ArrayList;
import java.util.List;

import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;

/**
 * Deterministic double-precision reference for the Bermudan bond call: finite differences for the CIR pricing PDE
 * (in time to go {@code tau})
 *
 * <pre>
 *   V_tau = theta (mu - r) V_r + 1/2 sigma^2 r V_rr - r V
 * </pre>
 *
 * between exercise dates, with {@code V = max(V, E_j(r))} at each exercise date and {@code V = E_m(r)} at the last.
 *
 * <ul>
 *   <li>Rate grid on {@code [0, rMax]} with {@code r0} as a grid node: either uniform ({@link Grid#UNIFORM}) or
 *       stretched by {@code r = r0 + w sinh(xi)} with uniform {@code xi} ({@link Grid#SINH}), which concentrates
 *       nodes around {@code r0} (and the exercise boundaries next to it) and coarsens towards {@code rMax}.</li>
 *   <li>Second-order three-point differences on the non-uniform grid. Drift: central where both off-diagonal
 *       coefficients stay non-negative (discrete maximum principle), otherwise first-order upwind (only near
 *       {@code r = 0}).</li>
 *   <li>{@code r = 0}: the diffusion vanishes and the drift {@code theta mu > 0} points inwards, so the PDE itself is
 *       used with a one-sided (upwind) difference; no artificial boundary value is imposed.</li>
 *   <li>{@code r = rMax}: linearity condition {@code V_rr = 0}.</li>
 *   <li>Time: Crank-Nicolson with Rannacher start-up (4 implicit Euler half steps) after every exercise date,
 *       where the payoff has a kink.</li>
 *   <li>Kink treatment ({@code smoothKinks}): the kink of {@code max(C, E)} (and of the final payoff) lies between
 *       two nodes. Sampling it at the nodes gives an O(dr^2) error whose constant depends on where the kink sits
 *       inside its cell, so the price oscillates under refinement. With smoothing, the node whose dual cell
 *       {@code [r_{i-1/2}, r_{i+1/2}]} contains the kink receives the cell average of {@code max(C, E)} instead
 *       (continuation value by cubic interpolation, exercise value exact, kink located by bisection); all other
 *       nodes keep point values. This is the projection/averaging remedy of Pooley, Vetzal and Forsyth (2003).</li>
 * </ul>
 * The reference solves the continuous-time model (no Euler bias, no Monte Carlo noise, exact exercise policy on
 * the grid) but has its own discretisation error, estimated by grid refinement and by comparison with the
 * independent quadrature reference {@link CirQuadratureBermudan}.
 */
public final class CirFiniteDifferenceBermudan {

    public enum Grid { UNIFORM, SINH }

    /**
     * Discretisation scheme.
     *
     * @param grid        node placement
     * @param stretch     width {@code w} of the sinh stretching (smaller = stronger concentration at {@code r0});
     *                    ignored for {@link Grid#UNIFORM}
     * @param smoothKinks cell averaging at the kink of the exercise payoff
     * @param rateMax     upper end of the rate grid (approximate: the grid is adjusted so that {@code r0} is a node)
     */
    public record Scheme(Grid grid, double stretch, boolean smoothKinks, double rateMax) {
        public Scheme withRateMax(double r) { return new Scheme(grid, stretch, smoothKinks, r); }
        public Scheme withSmoothing(boolean s) { return new Scheme(grid, stretch, s, rateMax); }
        public String label() {
            return (grid == Grid.UNIFORM ? "uniform" : "sinh(w=" + stretch + ")") + (smoothKinks ? "+smooth" : "");
        }
    }

    /**
     * @param exerciseBoundaries per exercise date the rates where {@code E_j = C_j} (last date: {@code B = K});
     *                           empty if there is none on the grid
     */
    public record Result(double price, int rateNodes, int timeStepsPerYear, double rMax, double[][] exerciseBoundaries) {
    }

    private final double theta;
    private final double mu;
    private final double sigma;
    private final double r0;
    private final CirZeroBondDouble zeroBond;

    public CirFiniteDifferenceBermudan(CirParameters p) {
        this.theta = p.meanReversionSpeedDouble();
        this.mu = p.longTermMeanDouble();
        this.sigma = p.volatilityDouble();
        this.r0 = p.initialRateDouble();
        this.zeroBond = new CirZeroBondDouble(p, ZeroBondForm.LOG1P_RATIO);
    }

    public double defaultRateMax() {
        double level = Math.max(mu, r0);
        double stationarySd = sigma * Math.sqrt(level / (2.0 * theta));
        return Math.max(0.5, 4.0 * level + 15.0 * stationarySd);
    }

    /** Uniform grid, point sampling of the kink: the scheme used for the first reported results. */
    public Scheme legacyScheme() {
        return new Scheme(Grid.UNIFORM, Double.NaN, false, defaultRateMax());
    }

    /** Default: sinh grid concentrated at {@code r0}, kink smoothing. */
    public Scheme defaultScheme() {
        return new Scheme(Grid.SINH, defaultStretch(), true, defaultRateMax());
    }

    public double defaultStretch() {
        return 0.5 * Math.max(mu, r0);
    }

    public Result price(BermudanCallOption option, int rateIntervals, int timeStepsPerYear) {
        return price(option, rateIntervals, timeStepsPerYear, defaultScheme());
    }

    public Result price(BermudanCallOption option, int rateIntervals, int timeStepsPerYear, Scheme scheme) {
        double[] r = grid(rateIntervals, scheme);
        int r0Index = indexOf(r, r0);
        int n = r.length - 1;

        // operator coefficients L V_i = a_i V_{i-1} + b_i V_i + c_i V_{i+1}
        double[] a = new double[n + 1];
        double[] b = new double[n + 1];
        double[] c = new double[n + 1];
        double s2 = sigma * sigma;
        b[0] = -theta * mu / (r[1] - r[0]);
        c[0] = theta * mu / (r[1] - r[0]);
        for (int i = 1; i < n; i++) {
            double hm = r[i] - r[i - 1];
            double hp = r[i + 1] - r[i];
            double drift = theta * (mu - r[i]);
            double d = 0.5 * s2 * r[i];
            double ai = 2.0 * d / (hm * (hm + hp));
            double ci = 2.0 * d / (hp * (hm + hp));
            double bi = -2.0 * d / (hm * hp) - r[i];
            double ac = -drift * hp / (hm * (hm + hp));
            double cc = drift * hm / (hp * (hm + hp));
            if (ai + ac >= 0.0 && ci + cc >= 0.0) {
                ai += ac;
                ci += cc;
                bi += drift * (hp - hm) / (hm * hp);
            } else if (drift > 0) {
                bi -= drift / hp;
                ci += drift / hp;
            } else {
                ai -= drift / hm;
                bi += drift / hm;
            }
            a[i] = ai;
            b[i] = bi;
            c[i] = ci;
        }
        // linear extrapolation V_n = V_{n-1} + q (V_{n-1} - V_{n-2})
        double q = (r[n] - r[n - 1]) / (r[n - 1] - r[n - 2]);

        RemainingBondDouble bond = new RemainingBondDouble(option, zeroBond);
        int m = option.exerciseDateCount();
        double[][] boundaries = new double[m][];
        double[] v = new double[n + 1];
        boundaries[m - 1] = terminalPayoff(r, v, bond, m - 1, scheme.smoothKinks());
        for (int j = m - 1; j >= 0; j--) {
            double from = option.exerciseDateDouble(j);
            double to = j == 0 ? 0.0 : option.exerciseDateDouble(j - 1);
            double length = from - to;
            int steps = Math.max(4, (int) Math.ceil(length * timeStepsPerYear));
            double dt = length / steps;
            // Rannacher: two full steps as four implicit half steps, then Crank-Nicolson
            for (int h = 0; h < 4; h++) {
                step(v, a, b, c, q, dt / 2.0, 1.0);
            }
            for (int s = 2; s < steps; s++) {
                step(v, a, b, c, q, dt, 0.5);
            }
            if (j > 0) {
                boundaries[j - 1] = exercise(r, v, bond, j - 1, scheme.smoothKinks());
            }
        }
        return new Result(v[r0Index], n, timeStepsPerYear, r[n], boundaries);
    }

    // ------------------------------------------------------------------------------------------------ grid

    double[] grid(int intervals, Scheme scheme) {
        if (intervals < 4) throw new IllegalArgumentException("need at least 4 rate intervals");
        double rMax = scheme.rateMax();
        if (scheme.grid() == Grid.UNIFORM) {
            int r0Index = Math.max(1, (int) Math.round(r0 / rMax * intervals));
            double dr = r0 / r0Index;
            double[] r = new double[intervals + 1];
            for (int i = 0; i <= intervals; i++) r[i] = i * dr;
            r[r0Index] = r0;
            return r;
        }
        double w = scheme.stretch();
        double lo = asinh(-r0 / w);
        double hi = asinh((rMax - r0) / w);
        int r0Index = Math.max(1, (int) Math.round(intervals * (-lo) / (hi - lo)));
        double dxi = -lo / r0Index;
        double[] r = new double[intervals + 1];
        for (int i = 0; i <= intervals; i++) {
            r[i] = r0 + w * Math.sinh(lo + i * dxi);
        }
        r[0] = 0.0;
        r[r0Index] = r0;
        return r;
    }

    private static double asinh(double x) {
        return Math.log(x + Math.sqrt(x * x + 1.0));
    }

    private static int indexOf(double[] r, double x) {
        for (int i = 0; i < r.length; i++) if (r[i] == x) return i;
        throw new IllegalStateException("r0 is not a grid node");
    }

    // ------------------------------------------------------------------------------------------------ payoff

    /** Final payoff {@code (B - K)^+}; returns the kink {@code B = K} if it lies inside the grid. */
    private double[] terminalPayoff(double[] r, double[] v, RemainingBondDouble bond, int j, boolean smooth) {
        int n = r.length - 1;
        for (int i = 0; i <= n; i++) {
            v[i] = bond.exerciseValue(bond.bondValue(j, r[i]));
        }
        List<Double> kinks = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            double gi = bond.bondValue(j, r[i]) - bond.strike();
            double gk = bond.bondValue(j, r[i + 1]) - bond.strike();
            if ((gi > 0) != (gk > 0)) {
                final int idx = i;
                double x = bisect(t -> bond.bondValue(j, t) - bond.strike(), r[idx], r[idx + 1]);
                kinks.add(x);
                if (smooth) {
                    int k = ownerOf(r, i, x);
                    v[k] = cellAverage(r, k, x, t -> bond.exerciseValue(bond.bondValue(j, t)));
                }
            }
        }
        return kinks.stream().mapToDouble(Double::doubleValue).toArray();
    }

    /** {@code V = max(C, E_j)} at an exercise date; returns the rates where {@code C = E_j}. */
    private double[] exercise(double[] r, double[] v, RemainingBondDouble bond, int j, boolean smooth) {
        int n = r.length - 1;
        double[] e = new double[n + 1];
        for (int i = 0; i <= n; i++) e[i] = bond.exerciseValue(bond.bondValue(j, r[i]));
        double[] cont = v.clone();
        List<Double> kinks = new ArrayList<>();
        double[] averaged = new double[n + 1];
        boolean[] isAveraged = new boolean[n + 1];
        for (int i = 0; i < n; i++) {
            double gi = cont[i] - e[i];
            double gk = cont[i + 1] - e[i + 1];
            if ((gi > 0) != (gk > 0)) {
                int lo = Math.max(0, Math.min(i - 1, n - 3));
                Fn cubic = t -> lagrange(r, cont, lo, t);
                Fn exerciseValue = t -> bond.exerciseValue(bond.bondValue(j, t));
                double x = bisect(t -> cubic.at(t) - exerciseValue.at(t), r[i], r[i + 1]);
                kinks.add(x);
                if (smooth) {
                    int k = ownerOf(r, i, x);
                    averaged[k] = cellAverage(r, k, x, t -> Math.max(cubic.at(t), exerciseValue.at(t)));
                    isAveraged[k] = true;
                }
            }
        }
        for (int i = 0; i <= n; i++) {
            v[i] = isAveraged[i] ? averaged[i] : Math.max(cont[i], e[i]);
        }
        return kinks.stream().mapToDouble(Double::doubleValue).toArray();
    }

    @FunctionalInterface
    private interface Fn {
        double at(double x);
    }

    /** The node among {@code i, i+1} whose dual cell contains {@code x}. */
    private static int ownerOf(double[] r, int i, double x) {
        return x < 0.5 * (r[i] + r[i + 1]) ? i : i + 1;
    }

    /** Average of {@code f} over the dual cell of node {@code k}, split at the kink {@code x} (Gauss-Legendre). */
    private static double cellAverage(double[] r, int k, double x, Fn f) {
        int n = r.length - 1;
        double lo = k == 0 ? r[0] : 0.5 * (r[k - 1] + r[k]);
        double hi = k == n ? r[n] : 0.5 * (r[k] + r[k + 1]);
        double kink = Math.min(hi, Math.max(lo, x));
        return (gauss(f, lo, kink) + gauss(f, kink, hi)) / (hi - lo);
    }

    private static final double[] GL_X = {
            -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
            0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    private static final double[] GL_W = {
            0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
            0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

    private static double gauss(Fn f, double a, double b) {
        if (b <= a) return 0.0;
        double m = 0.5 * (a + b);
        double hw = 0.5 * (b - a);
        double s = 0.0;
        for (int i = 0; i < GL_X.length; i++) s += GL_W[i] * f.at(m + hw * GL_X[i]);
        return s * hw;
    }

    /** Cubic Lagrange interpolation through nodes {@code lo .. lo+3}. */
    private static double lagrange(double[] r, double[] y, int lo, double x) {
        double s = 0.0;
        for (int i = lo; i < lo + 4; i++) {
            double l = 1.0;
            for (int k = lo; k < lo + 4; k++) {
                if (k != i) l *= (x - r[k]) / (r[i] - r[k]);
            }
            s += l * y[i];
        }
        return s;
    }

    private static double bisect(Fn g, double a, double b) {
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

    // ------------------------------------------------------------------------------------------------ time step

    /**
     * One theta-scheme step {@code (I - w dt L) V_new = (I + (1-w) dt L) V_old}, unknowns 0..n-1,
     * {@code V_n = V_{n-1} + q (V_{n-1} - V_{n-2})}.
     */
    private static void step(double[] v, double[] a, double[] b, double[] c, double q, double dt, double w) {
        int n = v.length - 1;
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++) {
            double lv = b[i] * v[i] + c[i] * v[i + 1] + (i > 0 ? a[i] * v[i - 1] : 0.0);
            rhs[i] = v[i] + (1.0 - w) * dt * lv;
        }
        double[] lower = new double[n];
        double[] diag = new double[n];
        double[] upper = new double[n];
        for (int i = 0; i < n; i++) {
            double ai = i > 0 ? a[i] : 0.0;
            double bi = b[i];
            double ci = c[i];
            if (i == n - 1) {
                bi += (1.0 + q) * ci;
                ai -= q * ci;
                ci = 0.0;
            }
            lower[i] = -w * dt * ai;
            diag[i] = 1.0 - w * dt * bi;
            upper[i] = -w * dt * ci;
        }
        // Thomas algorithm
        for (int i = 1; i < n; i++) {
            double f = lower[i] / diag[i - 1];
            diag[i] -= f * upper[i - 1];
            rhs[i] -= f * rhs[i - 1];
        }
        v[n - 1] = rhs[n - 1] / diag[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            v[i] = (rhs[i] - upper[i] * v[i + 1]) / diag[i];
        }
        v[n] = v[n - 1] + q * (v[n - 1] - v[n - 2]);
    }
}

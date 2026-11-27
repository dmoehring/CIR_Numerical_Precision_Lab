package lab.cir.analytic;

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
 *   <li>Uniform rate grid {@code r_i = i dr} on {@code [0, rMax]}, chosen so that {@code r0} is a grid node.</li>
 *   <li>Drift: central differences where the cell Peclet condition {@code |drift| dr <= sigma^2 r} holds,
 *       otherwise first-order upwind (only near {@code r = 0}).</li>
 *   <li>{@code r = 0}: the diffusion vanishes and the drift {@code theta mu > 0} points inwards, so the PDE itself is
 *       used with a one-sided (upwind) difference; no artificial boundary value is imposed.</li>
 *   <li>{@code r = rMax}: linearity condition {@code V_rr = 0}.</li>
 *   <li>Time: Crank-Nicolson with Rannacher start-up (4 implicit Euler half steps) after every exercise date,
 *       where the payoff has a kink.</li>
 * </ul>
 * The reference solves the continuous-time model (no Euler bias, no Monte Carlo noise, exact exercise policy on
 * the grid) but has its own discretisation error, estimated by grid refinement.
 */
public final class CirFiniteDifferenceBermudan {

    public record Result(double price, int rateNodes, int timeStepsPerYear, double rMax) {
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
        this.zeroBond = new CirZeroBondDouble(p);
    }

    public double defaultRateMax() {
        double level = Math.max(mu, r0);
        double stationarySd = sigma * Math.sqrt(level / (2.0 * theta));
        return Math.max(0.5, 4.0 * level + 15.0 * stationarySd);
    }

    public Result price(BermudanCallOption option, int rateIntervals, int timeStepsPerYear) {
        return price(option, rateIntervals, timeStepsPerYear, defaultRateMax());
    }

    public Result price(BermudanCallOption option, int rateIntervals, int timeStepsPerYear, double rateMaxHint) {
        int r0Index = Math.max(1, (int) Math.round(r0 / rateMaxHint * rateIntervals));
        double dr = r0 / r0Index;
        int n = rateIntervals;
        double rMax = n * dr;
        double[] r = new double[n + 1];
        for (int i = 0; i <= n; i++) r[i] = i * dr;

        // operator coefficients L V_i = a_i V_{i-1} + b_i V_i + c_i V_{i+1}
        double[] a = new double[n + 1];
        double[] b = new double[n + 1];
        double[] c = new double[n + 1];
        double s2 = sigma * sigma;
        b[0] = -theta * mu / dr;
        c[0] = theta * mu / dr;
        for (int i = 1; i <= n; i++) {
            double drift = theta * (mu - r[i]);
            double diff = 0.5 * s2 * r[i] / (dr * dr);
            a[i] = diff;
            b[i] = -2.0 * diff - r[i];
            c[i] = diff;
            if (Math.abs(drift) * dr <= s2 * r[i]) {
                a[i] -= drift / (2.0 * dr);
                c[i] += drift / (2.0 * dr);
            } else if (drift > 0) {
                b[i] -= drift / dr;
                c[i] += drift / dr;
            } else {
                a[i] += drift / dr;
                b[i] -= drift / dr;
            }
        }

        RemainingBondDouble bond = new RemainingBondDouble(option, zeroBond);
        int m = option.exerciseDateCount();
        double[] v = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            v[i] = bond.exerciseValue(bond.bondValue(m - 1, r[i]));
        }
        for (int j = m - 1; j >= 0; j--) {
            double from = option.exerciseDateDouble(j);
            double to = j == 0 ? 0.0 : option.exerciseDateDouble(j - 1);
            double length = from - to;
            int steps = Math.max(4, (int) Math.ceil(length * timeStepsPerYear));
            double dt = length / steps;
            // Rannacher: two full steps as four implicit half steps, then Crank-Nicolson
            for (int h = 0; h < 4; h++) {
                step(v, a, b, c, dt / 2.0, 1.0);
            }
            for (int s = 2; s < steps; s++) {
                step(v, a, b, c, dt, 0.5);
            }
            if (j > 0) {
                for (int i = 0; i <= n; i++) {
                    double e = bond.exerciseValue(bond.bondValue(j - 1, r[i]));
                    if (e > v[i]) v[i] = e;
                }
            }
        }
        return new Result(v[r0Index], n, timeStepsPerYear, rMax);
    }

    /** One theta-scheme step {@code (I - w dt L) V_new = (I + (1-w) dt L) V_old}, unknowns 0..n-1, V_n by linearity. */
    private static void step(double[] v, double[] a, double[] b, double[] c, double dt, double w) {
        int n = v.length - 1;
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++) {
            double lv = b[i] * v[i] + c[i] * v[i + 1] + (i > 0 ? a[i] * v[i - 1] : 0.0);
            rhs[i] = v[i] + (1.0 - w) * dt * lv;
        }
        // matrix rows; last unknown row uses V_n = 2 V_{n-1} - V_{n-2}
        double[] lower = new double[n];
        double[] diag = new double[n];
        double[] upper = new double[n];
        for (int i = 0; i < n; i++) {
            double ai = i > 0 ? a[i] : 0.0;
            double bi = b[i];
            double ci = c[i];
            if (i == n - 1) {
                bi += 2.0 * ci;
                ai -= ci;
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
        v[n] = 2.0 * v[n - 1] - v[n - 2];
    }
}

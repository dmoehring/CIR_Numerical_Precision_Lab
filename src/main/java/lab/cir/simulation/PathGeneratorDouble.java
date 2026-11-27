package lab.cir.simulation;

import lab.cir.analytic.RemainingBondDouble;
import lab.cir.model.CirParameters;
import lab.cir.random.NormalSource;

/**
 * Full-truncation Euler scheme for CIR in {@code double}:
 *
 * <pre>
 *   x_{k+1} = x_k + theta (mu - x_k^+) dt + sigma sqrt(x_k^+) dW_k,      x_0 = r_0
 *   I_{k+1} = I_k + x_k^+ dt                                            (left rectangle rule)
 * </pre>
 *
 * The internal state {@code x_k} may become negative; it is not projected. The effective rate {@code x_k^+}
 * is used for drift, diffusion, integration, bond valuation and clustering.
 */
public final class PathGeneratorDouble {

    private final double theta;
    private final double mu;
    private final double sigma;
    private final double r0;
    private final double dt;
    private final double sqrtStepDt;
    private final TimeGrid grid;

    public PathGeneratorDouble(CirParameters p, TimeGrid grid) {
        this.theta = p.meanReversionSpeedDouble();
        this.mu = p.longTermMeanDouble();
        this.sigma = p.volatilityDouble();
        this.r0 = p.initialRateDouble();
        this.grid = grid;
        this.dt = 1.0 / grid.stepsPerYear();
        // sqrt(dt) for m = 1, sqrt(dt / m) for coupled fine increments
        this.sqrtStepDt = grid.normalsPerStep() == 1 ? Math.sqrt(dt) : Math.sqrt(dt / grid.normalsPerStep());
    }

    public PathDataDouble simulate(NormalSource source, int paths, int threads) {
        checkSource(source);
        PathDataDouble data = new PathDataDouble(grid.exerciseDateCount(), paths);
        ParallelPaths.forEach(paths, threads, p -> {
            double[] z = source.path(p);
            data.negativeStateCount[p] = simulatePath(z, p, data, null, null);
        });
        return data;
    }

    public void computeExerciseValues(PathDataDouble data, RemainingBondDouble bond, int threads) {
        ParallelPaths.forEach(data.paths(), threads, p -> {
            for (int j = 0; j < data.dates(); j++) {
                double b = bond.bondValue(j, data.rate[j][p]);
                data.bondValue[j][p] = b;
                data.exerciseValue[j][p] = bond.exerciseValue(b);
            }
        });
    }

    /**
     * Full trajectory of one path (diagnostics only): effective rates {@code x_k^+} and integrals {@code I_k},
     * {@code k = 0..N}. Uses exactly the same code as {@link #simulate}.
     */
    public double[][] trace(NormalSource source, int pathIndex) {
        checkSource(source);
        double[] rates = new double[grid.numberOfSteps() + 1];
        double[] integrals = new double[grid.numberOfSteps() + 1];
        simulatePath(source.path(pathIndex), -1, null, rates, integrals);
        return new double[][]{rates, integrals};
    }

    private int simulatePath(double[] z, int p, PathDataDouble out, double[] traceRate, double[] traceIntegral) {
        final int n = grid.numberOfSteps();
        final int m = grid.normalsPerStep();
        double x = r0;
        double integral = 0.0;
        int j = 0;
        int negative = 0;
        for (int k = 0; k <= n; k++) {
            if (x < 0.0) {
                negative++;
            }
            double xp = x > 0.0 ? x : 0.0;
            if (traceRate != null) {
                traceRate[k] = xp;
                traceIntegral[k] = integral;
            }
            if (out != null && k == grid.exerciseStepIndex(j)) {
                out.rate[j][p] = xp;
                out.cumulativeIntegral[j][p] = integral;
                j++;
            }
            if (k == n) {
                break;
            }
            integral += xp * dt;
            double dW;
            if (m == 1) {
                dW = sqrtStepDt * z[k];
            } else {
                double s = 0.0;
                for (int i = 0; i < m; i++) {
                    s += z[k * m + i];
                }
                dW = sqrtStepDt * s;
            }
            x = x + theta * (mu - xp) * dt + sigma * Math.sqrt(xp) * dW;
        }
        return negative;
    }

    private void checkSource(NormalSource source) {
        if (source.normalsPerPath() != grid.normalsPerPath()) {
            throw new IllegalArgumentException("normal source does not match the time grid");
        }
    }

    public TimeGrid grid() { return grid; }
}

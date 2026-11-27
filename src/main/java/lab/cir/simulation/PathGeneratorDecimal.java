package lab.cir.simulation;

import lab.cir.analytic.RemainingBondDecimal;
import lab.cir.model.CirParameters;
import lab.cir.numeric.DecimalMath;
import lab.cir.random.NormalSource;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Full-truncation Euler scheme for CIR in {@code BigDecimal}; same scheme and order of operations as
 * {@link PathGeneratorDouble}.
 *
 * <p>Random input: the shared binary64 normals are converted <em>exactly</em> with {@code new BigDecimal(z)}.
 * The MathContext is first applied to random input in the multiplication {@code sqrt(dt) * Z}
 * (for {@code normalsPerStep > 1}: the exact decimal sum of the fine normals is multiplied by {@code sqrt(dt/m)}).
 */
public final class PathGeneratorDecimal {

    private final BigDecimal theta;
    private final BigDecimal mu;
    private final BigDecimal sigma;
    private final BigDecimal r0;
    private final BigDecimal dt;
    private final BigDecimal sqrtStepDt;
    private final TimeGrid grid;
    private final DecimalMath math;
    private final MathContext mc;

    public PathGeneratorDecimal(CirParameters p, TimeGrid grid, DecimalMath math) {
        this.math = math;
        this.mc = math.mc();
        this.theta = p.meanReversionSpeed();
        this.mu = p.longTermMean();
        this.sigma = p.volatility();
        this.r0 = p.initialRate();
        this.grid = grid;
        this.dt = BigDecimal.ONE.divide(BigDecimal.valueOf(grid.stepsPerYear()), mc);
        this.sqrtStepDt = grid.normalsPerStep() == 1 ? math.sqrt(dt)
                : math.sqrt(dt.divide(BigDecimal.valueOf(grid.normalsPerStep()), mc));
    }

    public PathDataDecimal simulate(NormalSource source, int paths, int threads) {
        checkSource(source);
        PathDataDecimal data = new PathDataDecimal(grid.exerciseDateCount(), paths);
        ParallelPaths.forEach(paths, threads, p -> {
            double[] z = source.path(p);
            data.negativeStateCount[p] = simulatePath(z, p, data, null, null);
        });
        return data;
    }

    public void computeExerciseValues(PathDataDecimal data, RemainingBondDecimal bond, int threads) {
        ParallelPaths.forEach(data.paths(), threads, p -> {
            for (int j = 0; j < data.dates(); j++) {
                BigDecimal b = bond.bondValue(j, data.rate[j][p]);
                data.bondValue[j][p] = b;
                data.exerciseValue[j][p] = bond.exerciseValue(b);
            }
        });
    }

    /** Full trajectory of one path for diagnostics (same code path as {@link #simulate}). */
    public BigDecimal[][] trace(NormalSource source, int pathIndex) {
        checkSource(source);
        BigDecimal[] rates = new BigDecimal[grid.numberOfSteps() + 1];
        BigDecimal[] integrals = new BigDecimal[grid.numberOfSteps() + 1];
        simulatePath(source.path(pathIndex), -1, null, rates, integrals);
        return new BigDecimal[][]{rates, integrals};
    }

    private int simulatePath(double[] z, int p, PathDataDecimal out, BigDecimal[] traceRate, BigDecimal[] traceIntegral) {
        final int n = grid.numberOfSteps();
        final int m = grid.normalsPerStep();
        BigDecimal x = r0;
        BigDecimal integral = BigDecimal.ZERO;
        int j = 0;
        int negative = 0;
        for (int k = 0; k <= n; k++) {
            if (x.signum() < 0) {
                negative++;
            }
            BigDecimal xp = DecimalMath.positivePart(x);
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
            integral = integral.add(xp.multiply(dt, mc), mc);
            BigDecimal dW;
            if (m == 1) {
                dW = sqrtStepDt.multiply(new BigDecimal(z[k]), mc);
            } else {
                BigDecimal s = BigDecimal.ZERO;
                for (int i = 0; i < m; i++) {
                    s = s.add(new BigDecimal(z[k * m + i])); // exact sum of exact binary values
                }
                dW = sqrtStepDt.multiply(s, mc);
            }
            BigDecimal drift = theta.multiply(mu.subtract(xp, mc), mc).multiply(dt, mc);
            BigDecimal diffusion = sigma.multiply(math.sqrt(xp), mc).multiply(dW, mc);
            x = x.add(drift, mc).add(diffusion, mc);
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

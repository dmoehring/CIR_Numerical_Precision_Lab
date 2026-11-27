package lab.cir.metrics;

import lab.cir.numeric.DecimalMath;
import lab.cir.random.NormalSource;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.simulation.PathGeneratorDecimal;
import lab.cir.simulation.PathGeneratorDouble;
import lab.cir.stats.DifferenceStats;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Differences of simulated quantities under common random numbers (12.2, 12.3):
 * at every exercise date over all paths (effective rate, cumulative integral, discount factor {@code exp(-I)},
 * bond value, exercise value), and along the full time grid for a small number of traced paths.
 * Discount factors are recomputed here for diagnostics; the valuation itself recomputes them on the fly.
 */
public final class StateDifferences {

    public record DateRow(int date, DifferenceStats rate, DifferenceStats integral, DifferenceStats discount,
                          DifferenceStats bond, DifferenceStats exercise) {
    }

    public record StepRow(int step, double time, DifferenceStats rate, DifferenceStats integral) {
    }

    public final DateRow[] dates;
    public final StepRow[] steps;

    public StateDifferences(PathDataDouble dd, PathDataDecimal dc, DecimalMath math) {
        this(dd, dc, math, null, null, null, 0);
    }

    public StateDifferences(PathDataDouble dd, PathDataDecimal dc, DecimalMath math, NormalSource source,
                            PathGeneratorDouble gd, PathGeneratorDecimal gc, int tracePaths) {
        MathContext mc = math.mc();
        dates = new DateRow[dd.dates()];
        for (int j = 0; j < dd.dates(); j++) {
            DifferenceStats r = new DifferenceStats(1e-12, mc);
            DifferenceStats in = new DifferenceStats(1e-12, mc);
            DifferenceStats df = new DifferenceStats(1e-12, mc);
            DifferenceStats b = new DifferenceStats(1e-12, mc);
            DifferenceStats e = new DifferenceStats(1e-9, mc);
            for (int p = 0; p < dd.paths(); p++) {
                r.add(dc.rate[j][p], dd.rate[j][p]);
                in.add(dc.cumulativeIntegral[j][p], dd.cumulativeIntegral[j][p]);
                BigDecimal dfc = math.exp(dc.cumulativeIntegral[j][p].negate());
                df.add(dfc, Math.exp(-dd.cumulativeIntegral[j][p]));
                b.add(dc.bondValue[j][p], dd.bondValue[j][p]);
                e.add(dc.exerciseValue[j][p], dd.exerciseValue[j][p]);
            }
            dates[j] = new DateRow(j, r, in, df, b, e);
        }
        if (source == null || tracePaths <= 0) {
            steps = new StepRow[0];
            return;
        }
        int n = gd.grid().numberOfSteps();
        steps = new StepRow[n + 1];
        for (int k = 0; k <= n; k++) {
            steps[k] = new StepRow(k, (double) k / gd.grid().stepsPerYear(), new DifferenceStats(1e-12, mc),
                    new DifferenceStats(1e-12, mc));
        }
        int traced = Math.min(tracePaths, dd.paths());
        for (int p = 0; p < traced; p++) {
            double[][] td = gd.trace(source, p);
            BigDecimal[][] tc = gc.trace(source, p);
            for (int k = 0; k <= n; k++) {
                steps[k].rate().add(tc[0][k], td[0][k]);
                steps[k].integral().add(tc[1][k], td[1][k]);
            }
        }
    }
}

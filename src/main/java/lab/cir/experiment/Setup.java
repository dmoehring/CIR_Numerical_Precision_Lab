package lab.cir.experiment;

import lab.cir.analytic.CirZeroBondDecimal;
import lab.cir.analytic.CirZeroBondDouble;
import lab.cir.analytic.RemainingBondDecimal;
import lab.cir.analytic.RemainingBondDouble;
import lab.cir.analytic.ZeroBondForm;
import lab.cir.config.ExperimentConfig;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;
import lab.cir.numeric.DecimalMath;
import lab.cir.numeric.PrecisionSettings;
import lab.cir.random.NormalSource;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.simulation.PathGeneratorDecimal;
import lab.cir.simulation.PathGeneratorDouble;
import lab.cir.simulation.TimeGrid;

/** All derived objects of one model/contract/grid/precision combination, for both arithmetics. */
public final class Setup {

    public final CirParameters model;
    public final BermudanCallOption option;
    public final TimeGrid grid;
    public final DecimalMath math;
    public final CirZeroBondDouble zeroBondDouble;
    public final CirZeroBondDecimal zeroBondDecimal;
    public final RemainingBondDouble bondDouble;
    public final RemainingBondDecimal bondDecimal;
    public final PathGeneratorDouble generatorDouble;
    public final PathGeneratorDecimal generatorDecimal;

    public Setup(CirParameters model, BermudanCallOption option, int stepsPerYear, int normalsPerStep, PrecisionSettings precision) {
        this(model, option, stepsPerYear, normalsPerStep, precision, ZeroBondForm.LOG_DIFFERENCE);
    }

    public Setup(CirParameters model, BermudanCallOption option, int stepsPerYear, int normalsPerStep, PrecisionSettings precision,
                 ZeroBondForm form) {
        this.model = model;
        this.option = option;
        this.grid = new TimeGrid(option, stepsPerYear, normalsPerStep);
        this.math = new DecimalMath(precision);
        this.zeroBondDouble = new CirZeroBondDouble(model, form);
        this.zeroBondDecimal = new CirZeroBondDecimal(model, math, form);
        this.bondDouble = new RemainingBondDouble(option, zeroBondDouble);
        this.bondDecimal = new RemainingBondDecimal(option, zeroBondDecimal);
        this.generatorDouble = new PathGeneratorDouble(model, grid);
        this.generatorDecimal = new PathGeneratorDecimal(model, grid, math);
    }

    public static Setup of(ExperimentConfig c) {
        return new Setup(c.model, c.option, c.stepsPerYear, c.normalsPerStep, c.precision, c.zeroBondForm);
    }

    public NormalSource source(long seed) {
        return new NormalSource(seed, grid.normalsPerPath());
    }

    public PathDataDouble simulateDouble(NormalSource src, int paths, int threads) {
        PathDataDouble d = generatorDouble.simulate(src, paths, threads);
        generatorDouble.computeExerciseValues(d, bondDouble, threads);
        return d;
    }

    public PathDataDecimal simulateDecimal(NormalSource src, int paths, int threads) {
        PathDataDecimal d = generatorDecimal.simulate(src, paths, threads);
        generatorDecimal.computeExerciseValues(d, bondDecimal, threads);
        return d;
    }
}

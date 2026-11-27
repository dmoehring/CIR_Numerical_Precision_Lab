package lab.cir.bench;

import lab.cir.analytic.CirZeroBondDecimal;
import lab.cir.analytic.CirZeroBondDouble;
import lab.cir.induction.BackwardInductionDecimal;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PathSelection;
import lab.cir.induction.RegressionVariable;
import lab.cir.induction.SummationMode;
import lab.cir.induction.TrainingStructure;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;
import lab.cir.model.CirParameters;
import lab.cir.numeric.PrecisionSettings;
import lab.cir.random.NormalSource;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.concurrent.TimeUnit;

/**
 * JMH benchmarks of the separate measurement areas (spec 18.1), double vs BigDecimal.
 * Inputs (paths, normals) are prepared in {@code @Setup} for the isolated areas; {@link #endToEndDouble} /
 * {@link #endToEndDecimal} include random number generation and conversions. Results are returned to JMH
 * (no dead-code elimination). Run with {@code -prof gc} for allocation rates.
 *
 * <pre>
 *   java -cp target/cir-numerical-precision-lab.jar org.openjdk.jmh.Main lab.cir.bench -prof gc
 * </pre>
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(value = 1, jvmArgsAppend = {"-Xmx2g"})
@State(Scope.Benchmark)
public class PrecisionBenchmarks {

    @Param({"1000"})
    public int paths;

    @Param({"16", "34", "64"})
    public int precision;

    private lab.cir.experiment.Setup setup;
    private NormalSource source;
    private PathDataDouble dataDouble;
    private PathDataDecimal dataDecimal;
    private TrainingStructure structureDouble;
    private TrainingStructure structureDecimal;
    private InductionSettings settings;
    private CirZeroBondDouble zbDouble;
    private CirZeroBondDecimal zbDecimal;
    private BigDecimal tau;
    private BigDecimal rate;

    @Setup(Level.Trial)
    public void setUp() {
        CirParameters model = CirParameters.of("0.50", "0.04", "0.10", "0.03");
        BermudanCallOption option = BermudanCallOption.of(CashflowSchedule.couponBond("100", "0.04", "5", 1), "100",
                "1", "2", "3", "4");
        setup = new lab.cir.experiment.Setup(model, option, 12, 1, new PrecisionSettings(precision, RoundingMode.HALF_EVEN, 10));
        source = setup.source(42);
        dataDouble = setup.simulateDouble(source, paths, 1);
        dataDecimal = setup.simulateDecimal(source, paths, 1);
        settings = new InductionSettings(ContinuationMethod.QUANTILE_CLUSTER, PathSelection.ALL, 10, 20,
                RegressionVariable.SHORT_RATE, SummationMode.NAIVE);
        structureDouble = TrainingStructure.fromDouble(dataDouble, settings);
        structureDecimal = TrainingStructure.fromDecimal(dataDecimal, settings);
        zbDouble = setup.zeroBondDouble;
        zbDecimal = setup.zeroBondDecimal;
        tau = new BigDecimal("3.0");
        rate = new BigDecimal("0.0312345");
    }

    // ---- analytic zero bond (coefficients + price)
    @Benchmark
    public double zeroBondDouble() {
        return zbDouble.price(3.0, 0.0312345);
    }

    @Benchmark
    public BigDecimal zeroBondDecimal() {
        return zbDecimal.price(tau, rate);
    }

    // ---- path simulation and pathwise integral (normals generated inside, as in production)
    @Benchmark
    public PathDataDouble simulationDouble() {
        return setup.generatorDouble.simulate(source, paths, 1);
    }

    @Benchmark
    public PathDataDecimal simulationDecimal() {
        return setup.generatorDecimal.simulate(source, paths, 1);
    }

    // ---- immediate exercise values (analytic bond values on all paths and dates)
    @Benchmark
    public PathDataDouble exerciseValuesDouble() {
        setup.generatorDouble.computeExerciseValues(dataDouble, setup.bondDouble, 1);
        return dataDouble;
    }

    @Benchmark
    public PathDataDecimal exerciseValuesDecimal() {
        setup.generatorDecimal.computeExerciseValues(dataDecimal, setup.bondDecimal, 1);
        return dataDecimal;
    }

    // ---- sorting and quantile partition
    @Benchmark
    public TrainingStructure clusteringDouble() {
        return TrainingStructure.fromDouble(dataDouble, settings);
    }

    @Benchmark
    public TrainingStructure clusteringDecimal() {
        return TrainingStructure.fromDecimal(dataDecimal, settings);
    }

    // ---- backward induction (targets, discounting, cluster means, decisions)
    @Benchmark
    public InductionResultDouble inductionDouble() {
        return BackwardInductionDouble.run(dataDouble, structureDouble, settings);
    }

    @Benchmark
    public InductionResultDecimal inductionDecimal() {
        return BackwardInductionDecimal.run(dataDecimal, structureDecimal, settings, setup.math);
    }

    // ---- end to end: normals, simulation, exercise values, clustering, induction
    @Benchmark
    public InductionResultDouble endToEndDouble() {
        PathDataDouble d = setup.simulateDouble(source, paths, 1);
        return BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, settings), settings);
    }

    @Benchmark
    public InductionResultDecimal endToEndDecimal() {
        PathDataDecimal d = setup.simulateDecimal(source, paths, 1);
        return BackwardInductionDecimal.run(d, TrainingStructure.fromDecimal(d, settings), settings, setup.math);
    }
}

package lab.cir.experiment;

import lab.cir.config.ExperimentConfig;
import lab.cir.induction.BackwardInductionDecimal;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PolicyValuation;
import lab.cir.induction.SummationMode;
import lab.cir.induction.TrainingStructure;
import lab.cir.induction.ValuationResultDecimal;
import lab.cir.induction.ValuationResultDouble;
import lab.cir.metrics.InductionDetails;
import lab.cir.metrics.PolicyComparison;
import lab.cir.metrics.PolicyView;
import lab.cir.metrics.StateDifferences;
import lab.cir.random.NormalSource;
import lab.cir.report.Measured;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.stats.PairedDifference;
import lab.cir.stats.SampleStatsDecimal;
import lab.cir.stats.SampleStatsDouble;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

/**
 * Main experiment: double vs BigDecimal valuation of the Bermudan option with common random numbers.
 *
 * <ol>
 *   <li>simulate training and (independent seed) valuation paths in both arithmetics;</li>
 *   <li>per continuation method: train the double policy; train the BigDecimal policy in Mode A (structure from
 *       double states) and Mode B (own structure); optional double+Kahan control;</li>
 *   <li>apply the fixed policies to the valuation paths (out-of-sample) and the frozen double policy with
 *       BigDecimal cashflows (pure arithmetic control);</li>
 *   <li>compute all comparison metrics.</li>
 * </ol>
 * Timings cover only the computation (no metrics, no I/O).
 */
public final class ComparisonExperiment {

    public static final class MethodRun {
        public InductionSettings settings;
        public Measured<TrainingStructure> structureDouble;
        public Measured<TrainingStructure> structureDecimal;
        public Measured<InductionResultDouble> inductionDouble;
        public Measured<InductionResultDecimal> inductionDecimalA;
        public Measured<InductionResultDecimal> inductionDecimalB;
        public InductionResultDouble inductionDoubleKahan;
        public Measured<ValuationResultDouble> valuationDouble;
        public Measured<ValuationResultDecimal> valuationDecimalA;
        public Measured<ValuationResultDecimal> valuationDecimalB;
        public ValuationResultDouble valuationDoubleKahan;
        /** Frozen double policy (out-of-sample exercise indices) with BigDecimal cashflows. */
        public PairedDifference frozenPolicy;
        public SampleStatsDecimal frozenPolicyPrice;
        public PolicyComparison inSampleA;
        public PolicyComparison inSampleB;
        public PolicyComparison outOfSampleA;
        public PolicyComparison outOfSampleB;
        public InductionDetails detailsA;
        public InductionDetails detailsB;

        public ContinuationMethod method() { return settings.method(); }
    }

    public static final class Run {
        public ExperimentConfig config;
        public Setup setup;
        public Measured<PathDataDouble> trainingDouble;
        public Measured<PathDataDecimal> trainingDecimal;
        public Measured<PathDataDouble> valuationDouble;
        public Measured<PathDataDecimal> valuationDecimal;
        public StateDifferences trainingStates;
        public final List<MethodRun> methods = new ArrayList<>();

        public MethodRun method(ContinuationMethod m) {
            return methods.stream().filter(r -> r.method() == m).findFirst().orElse(null);
        }
    }

    private ComparisonExperiment() {
    }

    public static Run run(ExperimentConfig cfg, boolean traceDiagnostics) {
        Run run = new Run();
        run.config = cfg;
        Setup s = Setup.of(cfg);
        run.setup = s;
        NormalSource training = s.source(cfg.trainingSeed);
        NormalSource valuation = s.source(cfg.valuationSeed);

        run.trainingDouble = Measured.of(() -> s.simulateDouble(training, cfg.trainingPaths, cfg.threads));
        run.trainingDecimal = Measured.of(() -> s.simulateDecimal(training, cfg.trainingPaths, cfg.threads));
        run.valuationDouble = Measured.of(() -> s.simulateDouble(valuation, cfg.valuationPaths, cfg.threads));
        run.valuationDecimal = Measured.of(() -> s.simulateDecimal(valuation, cfg.valuationPaths, cfg.threads));
        PathDataDouble td = run.trainingDouble.value();
        PathDataDecimal tc = run.trainingDecimal.value();
        PathDataDouble vd = run.valuationDouble.value();
        PathDataDecimal vc = run.valuationDecimal.value();

        for (ContinuationMethod method : cfg.methods) {
            MethodRun mr = new MethodRun();
            InductionSettings is = cfg.inductionSettings(method);
            mr.settings = is;
            mr.structureDouble = Measured.of(() -> TrainingStructure.fromDouble(td, is));
            mr.structureDecimal = Measured.of(() -> TrainingStructure.fromDecimal(tc, is));
            TrainingStructure sd = mr.structureDouble.value();
            TrainingStructure sc = mr.structureDecimal.value();
            mr.inductionDouble = Measured.of(() -> BackwardInductionDouble.run(td, sd, is));
            mr.inductionDecimalA = Measured.of(() -> BackwardInductionDecimal.run(tc, sd, is, s.math));
            mr.inductionDecimalB = Measured.of(() -> BackwardInductionDecimal.run(tc, sc, is, s.math));
            mr.valuationDouble = Measured.of(() -> PolicyValuation.valueDouble(vd, mr.inductionDouble.value()));
            mr.valuationDecimalA = Measured.of(() -> PolicyValuation.valueDecimal(vc, mr.inductionDecimalA.value(), vd, s.math));
            mr.valuationDecimalB = Measured.of(() -> PolicyValuation.valueDecimal(vc, mr.inductionDecimalB.value(), null, s.math));
            if (cfg.kahanControl && method == ContinuationMethod.QUANTILE_CLUSTER) {
                mr.inductionDoubleKahan = BackwardInductionDouble.run(td, sd, is.withSummation(SummationMode.KAHAN));
                mr.valuationDoubleKahan = PolicyValuation.valueDouble(vd, mr.inductionDoubleKahan);
            }
            ValuationResultDouble vdr = mr.valuationDouble.value();
            BigDecimal[] frozen = PolicyValuation.decimalCashflowsForFixedExercise(vc, vdr.exerciseIndex, s.math);
            mr.frozenPolicy = PairedDifference.of(frozen, vdr.pathValue, null, null, s.math);
            mr.frozenPolicyPrice = SampleStatsDecimal.of(frozen, s.math);

            mr.inSampleA = new PolicyComparison(PolicyView.OfDouble.of(mr.inductionDouble.value()),
                    PolicyView.OfDecimal.of(mr.inductionDecimalA.value()), td, tc, s.math);
            mr.inSampleB = new PolicyComparison(PolicyView.OfDouble.of(mr.inductionDouble.value()),
                    PolicyView.OfDecimal.of(mr.inductionDecimalB.value()), td, tc, s.math);
            mr.outOfSampleA = new PolicyComparison(PolicyView.OfDouble.of(vdr),
                    PolicyView.OfDecimal.of(mr.valuationDecimalA.value()), vd, vc, s.math);
            mr.outOfSampleB = new PolicyComparison(PolicyView.OfDouble.of(vdr),
                    PolicyView.OfDecimal.of(mr.valuationDecimalB.value()), vd, vc, s.math);
            mr.detailsA = new InductionDetails(mr.inductionDouble.value(), mr.inductionDecimalA.value(), td, tc, s.math);
            mr.detailsB = new InductionDetails(mr.inductionDouble.value(), mr.inductionDecimalB.value(), td, tc, s.math);
            run.methods.add(mr);
        }
        run.trainingStates = traceDiagnostics
                ? new StateDifferences(td, tc, s.math, training, s.generatorDouble, s.generatorDecimal, cfg.tracePaths)
                : new StateDifferences(td, tc, s.math);
        return run;
    }

    /** Fraction of paths with a (hypothetical) exercise decision at each date. */
    public static double[] exerciseRates(boolean[][] decision) {
        double[] r = new double[decision.length];
        for (int j = 0; j < decision.length; j++) {
            int c = 0;
            for (boolean b : decision[j]) if (b) c++;
            r[j] = (double) c / decision[j].length;
        }
        return r;
    }

    /** Fraction of paths with realised exercise at each date. */
    public static double[] realisedRates(int[] exerciseIndex, int dates) {
        double[] r = new double[dates];
        for (int t : exerciseIndex) if (t >= 0) r[t] += 1.0;
        for (int j = 0; j < dates; j++) r[j] /= exerciseIndex.length;
        return r;
    }

    static SampleStatsDouble stats(double[] v) {
        return SampleStatsDouble.of(v);
    }
}

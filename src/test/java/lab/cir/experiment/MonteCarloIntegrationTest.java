package lab.cir.experiment;

import lab.cir.TestFixtures;
import lab.cir.analytic.CirBondOptionAnalytic;
import lab.cir.config.ExperimentConfig;
import lab.cir.induction.BackwardInductionDecimal;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PolicyValuation;
import lab.cir.induction.TrainingStructure;
import lab.cir.induction.ValuationResultDouble;
import lab.cir.metrics.PolicyComparison;
import lab.cir.model.BermudanCallOption;
import lab.cir.numeric.PrecisionSettings;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.RoundingMode;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Monte Carlo properties with fixed seeds (deterministic: no sporadic failures). Statistical tolerances are stated
 * in units of the Monte Carlo standard error.
 */
class MonteCarloIntegrationTest {

    private static final InductionSettings CLUSTERS = InductionSettings.clusters(10, 20);

    private static InductionResultDouble inSample(Setup s, long seed, int paths) {
        PathDataDouble d = s.simulateDouble(s.source(seed), paths, 1);
        return BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, CLUSTERS), CLUSTERS);
    }

    @Test
    void singleExerciseDateIsTheEuropeanDiscountedPayoff() {
        BermudanCallOption eu = TestFixtures.bermudan().europeanAt(1);
        Setup s = new Setup(TestFixtures.model(), eu, 12, 1, PrecisionSettings.DECIMAL128);
        PathDataDouble d = s.simulateDouble(s.source(7), 2000, 1);
        InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, CLUSTERS), CLUSTERS);
        double sum = 0;
        for (int p = 0; p < 2000; p++) sum += d.exerciseValue[0][p] * Math.exp(-d.cumulativeIntegral[0][p]);
        assertEquals(sum / 2000, r.price.mean(), 1e-14);
    }

    @Test
    void europeanMonteCarloAgreesWithJamshidian() {
        CirBondOptionAnalytic a = new CirBondOptionAnalytic(TestFixtures.model());
        for (int j = 0; j < 4; j++) {
            BermudanCallOption eu = TestFixtures.bermudan().europeanAt(j);
            Setup s = new Setup(TestFixtures.model(), eu, 12, 1, PrecisionSettings.DECIMAL128);
            InductionResultDouble r = inSample(s, 4242, 20000);
            double ref = a.europeanCall(eu);
            double z = (r.price.mean() - ref) / r.price.standardError();
            assertTrue(Math.abs(z) < 4, "e" + (j + 1) + ": MC " + r.price.mean() + " ref " + ref + " z " + z);
        }
    }

    @Test
    void doubleAndDecimalAgreeClosely() {
        Setup s = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, PrecisionSettings.DECIMAL128);
        PathDataDouble d = s.simulateDouble(s.source(42), 300, 1);
        PathDataDecimal c = s.simulateDecimal(s.source(42), 300, 1);
        TrainingStructure sd = TrainingStructure.fromDouble(d, CLUSTERS);
        InductionResultDouble rd = BackwardInductionDouble.run(d, sd, CLUSTERS);
        InductionResultDecimal ra = BackwardInductionDecimal.run(c, sd, CLUSTERS, s.math);
        BigDecimal diff = ra.price.mean().subtract(new BigDecimal(rd.price.mean()));
        assertTrue(diff.abs().compareTo(new BigDecimal("1e-11")) < 0, "Delta V = " + diff);
        for (int j = 0; j < 4; j++) {
            assertArrayEquals(rd.cluster[j], ra.cluster[j], "Mode A: identical path indices per cluster");
        }
    }

    @Test
    void nonNegativeMonotoneInStrikeAndLinearScaling() {
        BermudanCallOption base = TestFixtures.bermudan().europeanAt(2);
        double prev = Double.POSITIVE_INFINITY;
        for (String k : new String[]{"96", "100", "104"}) {
            Setup s = new Setup(TestFixtures.model(), base.withStrike(new BigDecimal(k)), 12, 1, PrecisionSettings.DECIMAL128);
            InductionResultDouble r = inSample(s, 11, 1000);
            for (double v : r.pathValue) assertTrue(v >= 0.0);
            assertTrue(r.price.mean() <= prev, "European price non-increasing in K (pathwise)");
            prev = r.price.mean();
        }
        Setup s1 = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, PrecisionSettings.DECIMAL128);
        Setup s2 = new Setup(TestFixtures.model(), TestFixtures.bermudan().scaled(new BigDecimal("2.5")), 12, 1, PrecisionSettings.DECIMAL128);
        double v1 = inSample(s1, 11, 1000).price.mean();
        double v2 = inSample(s2, 11, 1000).price.mean();
        assertEquals(2.5 * v1, v2, 1e-12 * v2, "joint scaling of cashflows and strike");
    }

    @Test
    void bermudanOutOfSampleIsNotBelowEuropeanBeyondStatisticalTolerance() {
        Setup s = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, PrecisionSettings.DECIMAL128);
        PathDataDouble train = s.simulateDouble(s.source(1), 4000, 1);
        PathDataDouble val = s.simulateDouble(s.source(2), 4000, 1);
        InductionResultDouble ind = BackwardInductionDouble.run(train, TrainingStructure.fromDouble(train, CLUSTERS), CLUSTERS);
        ValuationResultDouble v = PolicyValuation.valueDouble(val, ind);
        CirBondOptionAnalytic a = new CirBondOptionAnalytic(TestFixtures.model());
        for (int j = 0; j < 4; j++) {
            double eu = a.europeanCall(TestFixtures.bermudan().europeanAt(j));
            assertTrue(v.price.mean() > eu - 4 * v.price.standardError(), "e" + (j + 1));
        }
    }

    @Test
    void policyOnTrainingPathsReproducesInSamplePrice() {
        Setup s = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, PrecisionSettings.DECIMAL128);
        PathDataDouble d = s.simulateDouble(s.source(3), 500, 1);
        for (ContinuationMethod m : ContinuationMethod.values()) {
            InductionSettings is = CLUSTERS.withMethod(m);
            InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, is), is);
            ValuationResultDouble v = PolicyValuation.valueDouble(d, r);
            assertEquals(r.price.mean(), v.price.mean(), 0.0, m.name());
        }
    }

    @Test
    void higherPrecisionReferenceIsStable() {
        PrecisionSettings p64 = new PrecisionSettings(64, RoundingMode.HALF_EVEN, 10);
        PrecisionSettings p100 = new PrecisionSettings(100, RoundingMode.HALF_EVEN, 10);
        Setup a = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, p64);
        Setup b = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, p100);
        PathDataDecimal da = a.simulateDecimal(a.source(5), 60, 1);
        PathDataDecimal db = b.simulateDecimal(b.source(5), 60, 1);
        InductionSettings one = InductionSettings.clusters(2, 20);
        BigDecimal va = BackwardInductionDecimal.run(da, TrainingStructure.fromDecimal(da, one), one, a.math).price.mean();
        BigDecimal vb = BackwardInductionDecimal.run(db, TrainingStructure.fromDecimal(db, one), one, b.math).price.mean();
        assertTrue(va.subtract(vb).abs().compareTo(new BigDecimal("1e-58")) < 0, va + " vs " + vb);
    }

    @Test
    void fullComparisonRunIsConsistent() {
        ExperimentConfig cfg = TestFixtures.config(200);
        ComparisonExperiment.Run run = ComparisonExperiment.run(cfg, true);
        assertEquals(2, run.methods.size());
        for (ComparisonExperiment.MethodRun m : run.methods) {
            for (PolicyComparison pc : new PolicyComparison[]{m.inSampleA, m.inSampleB, m.outOfSampleA, m.outOfSampleB}) {
                for (PolicyComparison.DateFigures f : pc.dates) {
                    assertEquals(200, f.agree() + f.differ());
                    assertEquals(f.differ(), f.doubleExDecimalCont() + f.decimalExDoubleCont());
                }
                assertTrue(pc.priceDifference.abs().compareTo(new BigDecimal("1e-10")) < 0);
            }
        }
        assertEquals(cfg.option.exerciseDateCount(), run.trainingStates.dates.length);
        assertTrue(run.trainingStates.steps.length > 0);
    }
}

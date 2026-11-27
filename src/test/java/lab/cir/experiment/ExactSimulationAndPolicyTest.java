package lab.cir.experiment;

import lab.cir.TestFixtures;
import lab.cir.analytic.CirBondOptionAnalytic;
import lab.cir.analytic.CirQuadratureBermudan;
import lab.cir.analytic.CirZeroBondDouble;
import lab.cir.analytic.ZeroBondForm;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PathSelection;
import lab.cir.induction.RegressionRuleDouble;
import lab.cir.induction.TrainingStructure;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;
import lab.cir.simulation.ExactPathGenerator;
import lab.cir.simulation.PathDataDouble;
import lab.cir.stats.SampleStatsDouble;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Exact simulation under the terminal forward measure and exact policy valuation by quadrature (fixed seeds,
 * statistical tolerances in units of the Monte Carlo standard error).
 */
class ExactSimulationAndPolicyTest {

    private static final CirQuadratureBermudan.Settings FAST = new CirQuadratureBermudan.Settings(1.0, 32, 12, 1e-11);

    @Test
    void exactPathsReproduceZeroBondsAndEuropeanCalls() {
        CirParameters p = TestFixtures.model();
        BermudanCallOption berm = TestFixtures.bermudan();
        ExactPathGenerator gen = new ExactPathGenerator(p, berm);
        PathDataDouble d = gen.simulate(2024, 50_000, 4);
        CirZeroBondDouble zb = new CirZeroBondDouble(p, ZeroBondForm.LOG1P_RATIO);
        CirBondOptionAnalytic a = new CirBondOptionAnalytic(p);
        int paths = d.paths();
        for (int j = 0; j < berm.exerciseDateCount(); j++) {
            double e = berm.exerciseDateDouble(j);
            // zero bond maturing at 5: worth P(e_j, 5; r_j) at e_j
            double[] zcb = new double[paths];
            double[] call = new double[paths];
            for (int i = 0; i < paths; i++) {
                double deflator = Math.exp(-d.cumulativeIntegral[j][i]); // P(0,S) / P(e_j,S)
                zcb[i] = zb.price(5.0 - e, d.rate[j][i]) * deflator;
                call[i] = d.exerciseValue[j][i] * deflator;
            }
            SampleStatsDouble z = SampleStatsDouble.of(zcb);
            assertEquals(zb.price(5.0, 0.03), z.mean(), 4 * z.standardError(), "zero bond via exact paths at e" + (j + 1));
            SampleStatsDouble c = SampleStatsDouble.of(call);
            double ref = a.europeanCall(berm.europeanAt(j));
            assertEquals(ref, c.mean(), 4 * c.standardError(), "European call via exact paths at e" + (j + 1));
            assertTrue(c.standardError() < 0.02);
        }
    }

    @Test
    void exactPathsAreIndependentOfThreadCount() {
        ExactPathGenerator gen = new ExactPathGenerator(TestFixtures.model(), TestFixtures.bermudan());
        PathDataDouble a = gen.simulate(7, 2000, 1);
        PathDataDouble b = gen.simulate(7, 2000, 4);
        for (int j = 0; j < a.dates(); j++) {
            org.junit.jupiter.api.Assertions.assertArrayEquals(a.rate[j], b.rate[j]);
        }
    }

    @Test
    void optimalRuleReproducesThePriceAndShiftedThresholdsLoseQuadratically() {
        BermudanCallOption berm = TestFixtures.bermudan();
        CirQuadratureBermudan q = new CirQuadratureBermudan(TestFixtures.model());
        CirQuadratureBermudan.Result opt = q.price(berm, FAST);
        double vStar = opt.price();
        int m = berm.exerciseDateCount();
        double[] rStar = new double[m - 1];
        for (int j = 0; j < m - 1; j++) rStar[j] = opt.exerciseBoundaries()[j][0];

        // "exercise iff E > C*" with the exact continuation function
        CirQuadratureBermudan.ExercisePolicy optimalRule =
                (j, r, b, e) -> e > 0 && (j == m - 1 || e > opt.continuationValue(j, r));
        assertEquals(vStar, q.valuePolicy(berm, optimalRule, FAST).price(), 1e-12);

        // threshold rule at the exact boundaries is optimal; shifted thresholds lose, quadratically in the shift
        double[] loss = new double[3];
        double[] shifts = {0.0, 0.001, 0.002};
        for (int k = 0; k < shifts.length; k++) {
            double shift = shifts[k];
            CirQuadratureBermudan.ExercisePolicy threshold =
                    (j, r, b, e) -> e > 0 && (j == m - 1 || r < rStar[j] + shift);
            loss[k] = vStar - q.valuePolicy(berm, threshold, FAST).price();
        }
        assertEquals(0.0, loss[0], 1e-11, "threshold at r* is optimal");
        assertTrue(loss[1] > 0 && loss[2] > loss[1], "shifted thresholds lose value");
        assertEquals(4.0, loss[2] / loss[1], 0.3, "loss ~ quadratic in the boundary shift");
    }

    @Test
    void trainedPoliciesNeverBeatTheOptimum() {
        BermudanCallOption berm = TestFixtures.bermudan();
        CirQuadratureBermudan q = new CirQuadratureBermudan(TestFixtures.model());
        double vStar = q.price(berm, FAST).price();
        PathDataDouble train = new ExactPathGenerator(TestFixtures.model(), berm).simulate(11, 2000, 1);
        InductionSettings cl = InductionSettings.clusters(10, 20);
        InductionSettings lsm = cl.withMethod(ContinuationMethod.LEAST_SQUARES).withPathSelection(PathSelection.IN_THE_MONEY)
                .withBasisDegree(3);
        for (InductionSettings s : new InductionSettings[]{cl, lsm}) {
            InductionResultDouble ind = BackwardInductionDouble.run(train, TrainingStructure.fromDouble(train, s), s);
            double v = q.valuePolicy(berm, MethodStudy.policyOf(ind), FAST).price();
            assertTrue(v <= vStar + 1e-12, s.label() + ": " + v + " > " + vStar);
            assertTrue(vStar - v < 0.05, s.label() + ": implausibly large loss " + (vStar - v));
        }
    }

    @Test
    void lsmBasisDegreeSetsTheNumberOfCoefficients() {
        PathDataDouble train = new ExactPathGenerator(TestFixtures.model(), TestFixtures.bermudan()).simulate(3, 500, 1);
        for (int degree = 1; degree <= 4; degree++) {
            InductionSettings s = InductionSettings.clusters(10, 20).withMethod(ContinuationMethod.LEAST_SQUARES)
                    .withBasisDegree(degree);
            InductionResultDouble ind = BackwardInductionDouble.run(train, TrainingStructure.fromDouble(train, s), s);
            assertEquals(degree + 1, ((RegressionRuleDouble) ind.rules[0]).coefficients().length);
        }
    }
}

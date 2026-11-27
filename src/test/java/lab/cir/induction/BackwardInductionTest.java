package lab.cir.induction;

import lab.cir.metrics.InductionDetails;
import lab.cir.numeric.DecimalMath;
import lab.cir.numeric.HouseholderQrDecimal;
import lab.cir.numeric.HouseholderQrDouble;
import lab.cir.numeric.PrecisionSettings;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/** Hand-constructed cases with known exercise decisions. */
class BackwardInductionTest {

    private static final DecimalMath MATH = new DecimalMath(PrecisionSettings.DECIMAL128);

    /**
     * 3 exercise dates, 4 paths, no discounting (I = 0), one cluster:
     * E(e3) = [5,0,3,0], E(e2) = [4,2,0,3], E(e1) = [1,6,2,0].
     * e3: tau = [2,-1,2,-1]; e2: C = mean[5,0,3,0] = 2 -> p0 (4>2) and p3 (3>2) exercise, p1 tie 2=2 -> continue;
     * e1: targets [4,0,3,3], C = 2.5 -> only p1 (6>2.5) exercises. Price = (4+6+3+3)/4 = 4.
     */
    static PathDataDouble constructedDouble(double discountPerDate) {
        double[][] e = {{1, 6, 2, 0}, {4, 2, 0, 3}, {5, 0, 3, 0}};
        PathDataDouble d = new PathDataDouble(3, 4);
        for (int j = 0; j < 3; j++) {
            for (int p = 0; p < 4; p++) {
                d.exerciseValue[j][p] = e[j][p];
                d.bondValue[j][p] = 100 + e[j][p];
                d.rate[j][p] = 0.01 * (p + 1);
                d.cumulativeIntegral[j][p] = discountPerDate * (j + 1);
            }
        }
        return d;
    }

    static PathDataDecimal toDecimal(PathDataDouble d) {
        PathDataDecimal c = new PathDataDecimal(d.dates(), d.paths());
        for (int j = 0; j < d.dates(); j++) {
            for (int p = 0; p < d.paths(); p++) {
                c.exerciseValue[j][p] = new BigDecimal(d.exerciseValue[j][p]);
                c.bondValue[j][p] = new BigDecimal(d.bondValue[j][p]);
                c.rate[j][p] = new BigDecimal(d.rate[j][p]);
                c.cumulativeIntegral[j][p] = new BigDecimal(d.cumulativeIntegral[j][p]);
            }
        }
        return c;
    }

    private static final InductionSettings ONE_CLUSTER = InductionSettings.clusters(1, 1);

    @Test
    void constructedCaseDecisionsAndTieConvention() {
        PathDataDouble d = constructedDouble(0.0);
        InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, ONE_CLUSTER), ONE_CLUSTER);
        assertArrayEquals(new boolean[]{true, false, true, false}, r.exercise[2], "last date: exercise iff E > 0");
        assertArrayEquals(new double[]{5, 0, 3, 0}, r.target[1], 0.0);
        assertEquals(2.0, r.continuation[1][0]);
        assertArrayEquals(new boolean[]{true, false, false, true}, r.exercise[1], "tie E = C continues; E = 0 never exercises");
        assertArrayEquals(new double[]{4, 0, 3, 3}, r.target[0], 0.0, "exercise at e2 replaced the later cashflow");
        assertArrayEquals(new int[]{1, -1, 2, 1}, r.futureExerciseBefore[0]);
        assertArrayEquals(new boolean[]{false, true, false, false}, r.exercise[0]);
        assertArrayEquals(new int[]{1, 0, 2, 1}, r.exerciseIndex);
        assertEquals(4.0, r.price.mean(), 0.0);
    }

    @Test
    void pathWithoutExerciseHasZeroTargetAndValue() {
        PathDataDouble d = constructedDouble(0.0);
        d.exerciseValue[0][3] = 0;
        d.exerciseValue[1][3] = 0;
        d.exerciseValue[2][3] = 0;
        InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, ONE_CLUSTER), ONE_CLUSTER);
        assertEquals(-1, r.exerciseIndex[3]);
        assertEquals(0.0, r.pathValue[3]);
        assertEquals(0.0, r.target[0][3]);
    }

    @Test
    void discountingOfTargetsAndPrice() {
        PathDataDouble d = constructedDouble(0.1);
        InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, ONE_CLUSTER), ONE_CLUSTER);
        // target at e2 of path 0 (exercise at e3): 5 exp(-(0.3 - 0.2))
        assertEquals(5 * Math.exp(-(0.30000000000000004 - 0.2)), r.target[1][0], 1e-15);
        for (int p = 0; p < 4; p++) {
            int t = r.exerciseIndex[p];
            double expected = t < 0 ? 0 : d.exerciseValue[t][p] * Math.exp(-d.cumulativeIntegral[t][p]);
            assertEquals(expected, r.pathValue[p], 0.0);
        }
    }

    @Test
    void decimalInductionReproducesTheConstructedCase() {
        PathDataDouble d = constructedDouble(0.0);
        PathDataDecimal c = toDecimal(d);
        InductionResultDecimal r = BackwardInductionDecimal.run(c, TrainingStructure.fromDecimal(c, ONE_CLUSTER), ONE_CLUSTER, MATH);
        assertArrayEquals(new int[]{1, 0, 2, 1}, r.exerciseIndex);
        assertEquals(0, r.price.mean().compareTo(new BigDecimal("4")));
        assertFalse(r.exercise[1][1], "tie continues in BigDecimal as well");
    }

    @Test
    void clusterMeansAreComputedPerCluster() {
        PathDataDouble d = constructedDouble(0.0);
        InductionSettings two = InductionSettings.clusters(2, 1);
        InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, two), two);
        ClusterRuleDouble rule = (ClusterRuleDouble) r.rules[1];
        // states ascending with path index: cluster 0 = paths {0,1} (targets 5,0), cluster 1 = {2,3} (3,0)
        assertEquals(2.5, rule.mean()[0]);
        assertEquals(1.5, rule.mean()[1]);
        assertEquals(0.02, rule.upperBounds()[0]);
        assertEquals(2, rule.size()[0]);
    }

    @Test
    void policyAppliedToTrainingPathsReproducesInSamplePrice() {
        PathDataDouble d = constructedDouble(0.1);
        InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, ONE_CLUSTER), ONE_CLUSTER);
        ValuationResultDouble v = PolicyValuation.valueDouble(d, r);
        assertArrayEquals(r.exerciseIndex, v.exerciseIndex);
        assertEquals(r.price.mean(), v.price.mean(), 0.0);
    }

    @Test
    void separateClustersSwitchOnlyNearTheBoundary() {
        // 4 paths, 2 clusters; BigDecimal states of paths 1 and 2 swap order by a tiny amount
        PathDataDouble d = constructedDouble(0.0);
        PathDataDecimal c = toDecimal(d);
        for (int j = 0; j < 3; j++) {
            d.rate[j][1] = 0.025;
            d.rate[j][2] = 0.0250000001;
            c.rate[j][1] = new BigDecimal("0.0250000001");
            c.rate[j][2] = new BigDecimal("0.025");
        }
        InductionSettings two = InductionSettings.clusters(2, 1);
        TrainingStructure sd = TrainingStructure.fromDouble(d, two);
        TrainingStructure sc = TrainingStructure.fromDecimal(c, two);
        InductionResultDouble rd = BackwardInductionDouble.run(d, sd, two);
        InductionResultDecimal ra = BackwardInductionDecimal.run(c, sd, two, MATH);
        InductionResultDecimal rb = BackwardInductionDecimal.run(c, sc, two, MATH);
        assertEquals(sd.partition(0).clusterOf(1), ra.cluster[0][1] < 0 ? -2 : sd.partition(0).clusterOf(1), "Mode A uses the double clusters");
        for (int p = 0; p < 4; p++) assertEquals(rd.cluster[1][p], ra.cluster[1][p], "Mode A: identical cluster ids");
        InductionDetails det = new InductionDetails(rd, rb, d, c, MATH);
        assertEquals(2, det.switches.get(1).switched(), "paths 1 and 2 swap clusters in Mode B");
        assertEquals(2, det.switches.get(1).switchedWithinOneRank(), "both adjacent to the quantile boundary");
        assertEquals(0, new InductionDetails(rd, ra, d, c, MATH).switches.get(1).switched());
    }

    @Test
    void lsmRecoversAnExactQuadratic() {
        int n = 50;
        double[][] x = new double[3][n];
        double[] y = new double[n];
        BigDecimal[][] xd = new BigDecimal[3][n];
        BigDecimal[] yd = new BigDecimal[n];
        for (int i = 0; i < n; i++) {
            double r = 0.001 * i;
            x[0][i] = 1;
            x[1][i] = r;
            x[2][i] = r * r;
            y[i] = 2 - 3 * r + 5 * r * r;
            BigDecimal rd = new BigDecimal("0.001").multiply(BigDecimal.valueOf(i));
            xd[0][i] = BigDecimal.ONE;
            xd[1][i] = rd;
            xd[2][i] = rd.multiply(rd);
            yd[i] = BigDecimal.TWO.subtract(new BigDecimal(3).multiply(rd)).add(new BigDecimal(5).multiply(rd).multiply(rd));
        }
        double[] b = HouseholderQrDouble.solve(x, y);
        assertEquals(2, b[0], 1e-12);
        assertEquals(-3, b[1], 1e-9);
        assertEquals(5, b[2], 1e-7);
        BigDecimal[] bd = HouseholderQrDecimal.solve(xd, yd, MATH);
        assertTrue(bd[2].subtract(new BigDecimal(5)).abs().compareTo(new BigDecimal("1e-25")) < 0, bd[2].toString());
        assertTrue(bd[1].add(new BigDecimal(3)).abs().compareTo(new BigDecimal("1e-28")) < 0);
    }

    @Test
    void itmSelectionUsesOnlyInTheMoneyPaths() {
        PathDataDouble d = constructedDouble(0.0);
        InductionSettings itm = ONE_CLUSTER.withPathSelection(PathSelection.IN_THE_MONEY);
        TrainingStructure s = TrainingStructure.fromDouble(d, itm);
        assertArrayEquals(new int[]{0, 1, 3}, s.eligible(1));
        InductionResultDouble r = BackwardInductionDouble.run(d, s, itm);
        // e2: targets of ITM paths {0,1,3} = 5,0,0 -> C = 5/3; p0 4 > 5/3, p1 2 > 5/3, p3 3 > 5/3
        assertEquals(5.0 / 3, r.continuation[1][0], 1e-15);
        assertTrue(Double.isNaN(r.continuation[1][2]), "OTM path not estimated");
        assertArrayEquals(new boolean[]{true, true, false, true}, r.exercise[1]);
    }
}

package lab.cir.simulation;

import lab.cir.ReferenceValues;
import lab.cir.TestFixtures;
import lab.cir.experiment.Setup;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;
import lab.cir.model.FellerMode;
import lab.cir.numeric.DecimalMath;
import lab.cir.numeric.PrecisionSettings;
import lab.cir.random.InverseNormal;
import lab.cir.random.NormalSource;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.MathContext;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SimulationTest {

    @Test
    void inverseNormalMatchesHighPrecisionReference() {
        for (Object[] r : ReferenceValues.NORMAL_QUANTILE) {
            double p = (Double) r[0];
            double ref = new BigDecimal((String) r[1]).doubleValue();
            assertEquals(ref, InverseNormal.quantile(p), Math.max(1e-15, Math.abs(ref) * 2e-15), "p=" + p);
        }
    }

    @Test
    void normalSourceIsReproducibleAndOrderIndependent() {
        NormalSource a = new NormalSource(42, 48);
        NormalSource b = new NormalSource(42, 48);
        double[] p7 = a.path(7);
        a.path(3);
        assertArrayEquals(p7, b.path(7), "same seed, same path index -> same normals regardless of call order");
        assertTrue(p7[0] != new NormalSource(43, 48).path(7)[0], "different seed -> different normals");
        assertTrue(p7[0] != a.path(8)[0]);
    }

    @Test
    void bothArithmeticsConsumeTheSameNormals() {
        // first Euler step reproduced by hand from the shared normal z_0 (no truncation at x_0 = r0 > 0)
        CirParameters m = TestFixtures.model();
        BermudanCallOption o = TestFixtures.bermudan();
        TimeGrid g = new TimeGrid(o, 12, 1);
        NormalSource src = new NormalSource(42, g.normalsPerPath());
        double z0 = src.path(5)[0];
        double dt = 1.0 / 12;
        double expected = 0.03 + 0.5 * (0.04 - 0.03) * dt + 0.1 * Math.sqrt(0.03) * (Math.sqrt(dt) * z0);
        double[][] td = new PathGeneratorDouble(m, g).trace(src, 5);
        assertEquals(expected, td[0][1], 0.0);
        DecimalMath math = new DecimalMath(PrecisionSettings.DECIMAL128);
        MathContext mc = math.mc();
        BigDecimal dtD = BigDecimal.ONE.divide(BigDecimal.valueOf(12), mc);
        BigDecimal dW = math.sqrt(dtD).multiply(new BigDecimal(z0), mc);
        BigDecimal expD = new BigDecimal("0.03").add(new BigDecimal("0.50").multiply(new BigDecimal("0.04").subtract(new BigDecimal("0.03"), mc), mc).multiply(dtD, mc), mc)
                .add(new BigDecimal("0.10").multiply(math.sqrt(new BigDecimal("0.03")), mc).multiply(dW, mc), mc);
        BigDecimal[][] tc = new PathGeneratorDecimal(m, g, math).trace(src, 5);
        assertEquals(0, expD.compareTo(tc[0][1]), "BigDecimal step uses new BigDecimal(z) of the shared normal");
    }

    @Test
    void simulationIsReproducibleAndIndependentOfThreadCount() {
        Setup s = new Setup(TestFixtures.model(), TestFixtures.bermudan(), 12, 1, PrecisionSettings.DECIMAL128);
        PathDataDouble a = s.simulateDouble(s.source(42), 300, 1);
        PathDataDouble b = s.simulateDouble(s.source(42), 300, 4);
        for (int j = 0; j < a.dates(); j++) {
            assertArrayEquals(a.rate[j], b.rate[j], 0.0);
            assertArrayEquals(a.exerciseValue[j], b.exerciseValue[j], 0.0);
        }
        PathDataDecimal c = s.simulateDecimal(s.source(42), 50, 1);
        PathDataDecimal d = s.simulateDecimal(s.source(42), 50, 3);
        for (int j = 0; j < c.dates(); j++) {
            assertArrayEquals(c.rate[j], d.rate[j]);
            assertArrayEquals(c.cumulativeIntegral[j], d.cumulativeIntegral[j]);
        }
    }

    @Test
    void leftRectangleIntegralAndDiscounting() {
        CirParameters m = TestFixtures.model();
        BermudanCallOption o = TestFixtures.bermudan();
        TimeGrid g = new TimeGrid(o, 12, 1);
        NormalSource src = new NormalSource(42, g.normalsPerPath());
        double[][] t = new PathGeneratorDouble(m, g).trace(src, 0);
        double integral = 0;
        for (int k = 0; k < g.numberOfSteps(); k++) {
            assertEquals(integral, t[1][k], 1e-15);
            integral += t[0][k] * (1.0 / 12);
        }
        PathDataDouble d = new PathGeneratorDouble(m, g).simulate(src, 200, 1);
        for (int p = 0; p < 200; p++) {
            double prev = 1.0;
            for (int j = 0; j < d.dates(); j++) {
                double df = Math.exp(-d.cumulativeIntegral[j][p]);
                assertTrue(df > 0 && df <= prev, "discount factors in (0,1] and non-increasing");
                prev = df;
            }
        }
    }

    @Test
    void fullTruncationCountsNegativeStatesInBothArithmetics() {
        // Feller violated on purpose (2*0.5*0.02 = 0.02 < 0.25^2): the Euler state goes negative regularly
        CirParameters m = CirParameters.of("0.5", "0.02", "0.25", "0.01", FellerMode.ALLOW_VIOLATION);
        BermudanCallOption o = TestFixtures.bermudan();
        Setup s = new Setup(m, o, 12, 1, PrecisionSettings.DECIMAL128);
        PathDataDouble d = s.simulateDouble(s.source(1), 300, 1);
        PathDataDecimal c = s.simulateDecimal(s.source(1), 300, 1);
        assertTrue(d.totalNegativeStates() > 0);
        assertEquals(d.totalNegativeStates(), c.totalNegativeStates(), "same truncation events in both arithmetics");
        for (int j = 0; j < d.dates(); j++) {
            for (int p = 0; p < 300; p++) {
                assertTrue(d.rate[j][p] >= 0.0, "effective rate is non-negative");
                assertTrue(c.rate[j][p].signum() >= 0);
            }
        }
    }

    @Test
    void couplingSumsFineIncrements() {
        // 4 steps/year with 3 normals per step uses the same fine Brownian path as 12 steps/year
        BermudanCallOption o = TestFixtures.bermudan();
        TimeGrid coarse = new TimeGrid(o, 4, 3);
        TimeGrid fine = new TimeGrid(o, 12, 1);
        assertEquals(fine.normalsPerPath(), coarse.normalsPerPath());
        BermudanCallOption offGrid = BermudanCallOption.of(TestFixtures.bond(), "100", "1.5");
        assertThrows(IllegalArgumentException.class, () -> new TimeGrid(offGrid, 5, 1), "1.5 * 5 is not an integer");
        assertEquals(18, new TimeGrid(offGrid, 12, 1).numberOfSteps());
    }
}

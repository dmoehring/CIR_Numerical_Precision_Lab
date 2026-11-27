package lab.cir.analytic;

import lab.cir.ReferenceValues;
import lab.cir.TestFixtures;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CirParameters;
import lab.cir.numeric.DecimalMath;
import lab.cir.numeric.PrecisionSettings;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class AnalyticTest {

    private static CirParameters params(Object[] row) {
        return CirParameters.of(BigDecimal.valueOf((Double) row[0]), BigDecimal.valueOf((Double) row[1]),
                BigDecimal.valueOf((Double) row[2]), new BigDecimal("0.03"), lab.cir.model.FellerMode.ALLOW_VIOLATION);
    }

    @Test
    void zeroBondDoubleMatchesIndependentReference() {
        for (ZeroBondForm form : ZeroBondForm.values()) {
            for (Object[] row : ReferenceValues.ZERO_BOND) {
                CirZeroBondDouble zb = new CirZeroBondDouble(params(row), form);
                double ref = new BigDecimal((String) row[5]).doubleValue();
                double p = zb.price((Double) row[3], (Double) row[4]);
                // error model: the absolute error of the log bracket (terms of size ~ h tau) is amplified by the exponent
                // c = 2 theta mu / sigma^2, so the relative error of P grows like eps * c * tau (measured: 3.8e-12 for c = 128,
                // tau = 30; this is a real limitation of the binary64 closed form, not of the test)
                double c = 2 * (Double) row[0] * (Double) row[1] / ((Double) row[2] * (Double) row[2]);
                double tol = 2e-15 * (1 + c * Math.max(1.0, (Double) row[3]));
                assertEquals(ref, p, Math.abs(ref) * tol, form + " " + java.util.Arrays.toString(row));
            }
        }
    }

    @Test
    void zeroBondDecimalMatchesIndependentReferenceTo30Digits() {
        DecimalMath math = new DecimalMath(PrecisionSettings.DECIMAL128);
        for (ZeroBondForm form : ZeroBondForm.values()) {
            for (Object[] row : ReferenceValues.ZERO_BOND) {
                // decimal inputs exactly as printed by the reference script
                CirParameters p = CirParameters.of(row[0].toString(), row[1].toString(), row[2].toString(), "0.03",
                        lab.cir.model.FellerMode.ALLOW_VIOLATION);
                CirZeroBondDecimal zb = new CirZeroBondDecimal(p, math, form);
                BigDecimal price = zb.price(new BigDecimal(row[3].toString()), new BigDecimal(row[4].toString()));
                BigDecimal ref = new BigDecimal((String) row[5]);
                BigDecimal rel = price.subtract(ref).abs().divide(ref, MathContext.DECIMAL64);
                double c = 2 * (Double) row[0] * (Double) row[1] / ((Double) row[2] * (Double) row[2]);
                // same error model at 34 digits: ~1e-33 * c * tau (exponent 40000 costs ~5 digits)
                double tol = 1e-32 * (1 + c * Math.max(1.0, (Double) row[3]));
                assertTrue(rel.doubleValue() < tol, form + " " + row[3] + ": rel " + rel);
            }
        }
    }

    @Test
    void zeroBondHigherPrecisionIsStable() {
        CirParameters p = TestFixtures.model();
        BigDecimal p64 = new CirZeroBondDecimal(p, new DecimalMath(new PrecisionSettings(64, RoundingMode.HALF_EVEN, 10)))
                .price(new BigDecimal("4"), new BigDecimal("0.03"));
        BigDecimal p100 = new CirZeroBondDecimal(p, new DecimalMath(new PrecisionSettings(100, RoundingMode.HALF_EVEN, 10)))
                .price(new BigDecimal("4"), new BigDecimal("0.03"));
        assertTrue(p64.subtract(p100).abs().compareTo(new BigDecimal("1e-61")) < 0, p64 + " vs " + p100);
    }

    @Test
    void smallVolatilityLimitIsDeterministicDiscounting() {
        // sigma -> 0: r(t) = mu + (r0 - mu) e^{-theta t}; P = exp(-mu tau - (r0 - mu)(1 - e^{-theta tau}) / theta); error O(sigma^2)
        CirParameters p = CirParameters.of("0.5", "0.04", "0.0001", "0.03");
        double tau = 5.0;
        double det = Math.exp(-(0.04 * tau + (0.03 - 0.04) * (1 - Math.exp(-0.5 * tau)) / 0.5));
        assertEquals(det, new CirZeroBondDouble(p).price(tau, 0.03), 1e-8);
    }

    @Test
    void zeroBondBasicProperties() {
        CirZeroBondDouble zb = new CirZeroBondDouble(TestFixtures.model());
        assertEquals(1.0, zb.price(0.0, 0.05));
        double prev = 1.0;
        for (double tau = 0.25; tau <= 30; tau += 0.25) {
            double v = zb.price(tau, 0.03);
            assertTrue(v < prev && v > 0);
            prev = v;
        }
        assertTrue(zb.price(2, 0.01) > zb.price(2, 0.05), "decreasing in r");
    }

    @Test
    void noncentralChiSquareMatchesReference() {
        for (double[] r : ReferenceValues.NONCENTRAL_CHI2) {
            double v = NoncentralChiSquared.cdf(r[0], r[1], r[2]);
            // validated range: relative error < 1e-12 for noncentrality <= 1000; for lambda ~ 1e4 the regularized gamma
            // function of Commons Math (a ~ 6000) limits the accuracy to ~1e-11 (measured 7e-12)
            double tol = r[2] <= 1000 ? 1e-12 : 2e-11;
            assertEquals(r[3], v, tol * r[3], java.util.Arrays.toString(r));
        }
        assertEquals(0.0, NoncentralChiSquared.cdf(0.0, 2, 1));
    }

    @Test
    void zeroBondCallMatchesReference() {
        for (double[] r : ReferenceValues.ZERO_BOND_CALL) {
            CirParameters p = CirParameters.of(BigDecimal.valueOf(r[0]), BigDecimal.valueOf(r[1]), BigDecimal.valueOf(r[2]),
                    BigDecimal.valueOf(r[3]), lab.cir.model.FellerMode.ALLOW_VIOLATION);
            double v = new CirBondOptionAnalytic(p).zeroBondCall(r[4], r[5], r[6]);
            assertEquals(r[7], v, 1e-12, java.util.Arrays.toString(r));
        }
    }

    @Test
    void couponBondCallJamshidianMatchesReference() {
        CirBondOptionAnalytic a = new CirBondOptionAnalytic(TestFixtures.model());
        for (double[] r : ReferenceValues.COUPON_BOND_CALL) {
            double v = a.couponBondCall(TestFixtures.bond(), r[0], r[1]);
            assertEquals(r[2], v, 1e-10, java.util.Arrays.toString(r));
            if (!Double.isNaN(r[3])) {
                assertEquals(r[3], a.criticalRate(TestFixtures.bond(), r[0], r[1]), 1e-12);
            } else {
                assertTrue(Double.isNaN(a.criticalRate(TestFixtures.bond(), r[0], r[1])), "no root: K >= bond value at r=0");
            }
        }
    }

    @Test
    void finiteDifferenceEuropeanMatchesJamshidianAndBermudanDominates() {
        CirParameters p = TestFixtures.model();
        BermudanCallOption berm = TestFixtures.bermudan();
        CirFiniteDifferenceBermudan fd = new CirFiniteDifferenceBermudan(p);
        CirBondOptionAnalytic a = new CirBondOptionAnalytic(p);
        double bermudan = fd.price(berm, 1000, 400).price();
        for (int j = 0; j < berm.exerciseDateCount(); j++) {
            double eu = fd.price(berm.europeanAt(j), 1000, 400).price();
            double ref = a.europeanCall(berm.europeanAt(j));
            assertEquals(ref, eu, 2e-4, "FD European vs Jamshidian at e" + (j + 1));
            assertTrue(bermudan >= eu - 1e-9, "Bermudan >= European at an exercise date");
        }
    }
}

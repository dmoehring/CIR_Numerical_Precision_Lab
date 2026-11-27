package lab.cir.model;

import lab.cir.TestFixtures;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ModelTest {

    @Test
    void fellerViolationIsRejectedInStrictMode() {
        // 2 * 0.5 * 0.04 = 0.04 < 0.3^2 = 0.09
        assertThrows(FellerConditionViolationException.class, () -> CirParameters.of("0.5", "0.04", "0.3", "0.03"));
    }

    @Test
    void fellerViolationIsAllowedOnlyExplicitly() {
        CirParameters p = CirParameters.of("0.5", "0.04", "0.3", "0.03", FellerMode.ALLOW_VIOLATION);
        assertFalse(p.satisfiesFellerCondition());
    }

    @Test
    void fellerBoundaryIsCheckedExactlyInDecimal() {
        // 2 * 0.5 * 0.04 = 0.04 = 0.2^2 exactly: accepted (>=), although 0.2*0.2 > 0.04 in binary64
        assertTrue(0.2 * 0.2 > 2 * 0.5 * 0.04, "binary64 would decide differently");
        assertDoesNotThrow(() -> CirParameters.of("0.5", "0.04", "0.2", "0.03"));
    }

    @Test
    void otherRestrictionsApplyInBothModes() {
        for (FellerMode mode : FellerMode.values()) {
            assertThrows(IllegalArgumentException.class, () -> CirParameters.of("0", "0.04", "0.1", "0.03", mode));
            assertThrows(IllegalArgumentException.class, () -> CirParameters.of("0.5", "-0.04", "0.1", "0.03", mode));
            assertThrows(IllegalArgumentException.class, () -> CirParameters.of("0.5", "0.04", "0", "0.03", mode));
            assertThrows(IllegalArgumentException.class, () -> CirParameters.of("0.5", "0.04", "0.1", "0", mode));
        }
    }

    @Test
    void doubleViewIsCorrectlyRoundedDecimalInput() {
        CirParameters p = TestFixtures.model();
        assertEquals(0.1, p.volatilityDouble());
        assertEquals(new BigDecimal("0.10"), p.volatility());
    }

    @Test
    void couponBondCashflows() {
        CashflowSchedule cf = TestFixtures.bond();
        assertEquals(5, cf.size());
        for (int i = 0; i < 4; i++) {
            assertEquals(0, cf.amount(i).compareTo(new BigDecimal("4")), "coupon c*N");
            assertEquals(0, cf.time(i).compareTo(BigDecimal.valueOf(i + 1)));
        }
        assertEquals(0, cf.amount(4).compareTo(new BigDecimal("104")), "(1+c) N at maturity");
        assertEquals(0, cf.time(4).compareTo(new BigDecimal("5")));
    }

    @Test
    void semiAnnualAndZeroCouponBonds() {
        CashflowSchedule semi = CashflowSchedule.couponBond("100", "0.05", "1", 2);
        assertEquals(2, semi.size());
        assertEquals(0, semi.amount(0).compareTo(new BigDecimal("2.5")));
        assertEquals(0, semi.time(0).compareTo(new BigDecimal("0.5")));
        CashflowSchedule zero = CashflowSchedule.couponBond("100", "0", "3", 1);
        assertEquals(1, zero.size(), "zero coupons are not positive cashflows");
        assertEquals(0, zero.amount(0).compareTo(new BigDecimal("100")));
    }

    @Test
    void cashflowAtExerciseDateIsNotPartOfRemainingBond() {
        CashflowSchedule cf = TestFixtures.bond();
        assertEquals(1, cf.firstIndexAfter(new BigDecimal("1.0")), "coupon at t=1 is paid, not bought");
        assertEquals(0, cf.firstIndexAfter(new BigDecimal("0.5")));
        assertEquals(4, cf.firstIndexAfter(new BigDecimal("4")), "only the final payment remains after t=4");
        assertEquals(0, cf.firstIndexAfter(BigDecimal.ZERO));
    }

    @Test
    void optionValidation() {
        CashflowSchedule cf = TestFixtures.bond();
        assertThrows(IllegalArgumentException.class, () -> BermudanCallOption.of(cf, "100", "2.0", "1.0"));
        assertThrows(IllegalArgumentException.class, () -> BermudanCallOption.of(cf, "100", "5.0"));
        assertThrows(IllegalArgumentException.class, () -> BermudanCallOption.of(cf, "0", "1.0"));
        assertThrows(IllegalArgumentException.class, () -> new CashflowSchedule(List.of(BigDecimal.ONE), List.of(BigDecimal.ZERO)));
        BermudanCallOption single = BermudanCallOption.of(CashflowSchedule.single("2", "1"), "0.9", "1");
        assertEquals(1, single.exerciseDateCount());
    }
}

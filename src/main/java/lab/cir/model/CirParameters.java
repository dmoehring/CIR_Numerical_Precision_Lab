package lab.cir.model;

import java.math.BigDecimal;
import java.util.Objects;

/**
 * Parameters of the CIR short-rate model under Q:
 * {@code dr = theta (mu - r) dt + sigma sqrt(r) dW}.
 *
 * <p>The canonical representation is decimal (parsed from strings). The {@code double} view is
 * derived from the same decimal string via {@link Double#parseDouble(String)}, i.e. it is the
 * correctly rounded binary64 value of the decimal input. Both numerical variants are therefore
 * derived from the same input configuration.
 *
 * <p>The Feller condition is checked exactly in decimal arithmetic (no rounding), so the
 * accept/reject decision does not depend on binary floating point.
 */
public final class CirParameters {

    private final BigDecimal meanReversionSpeed;
    private final BigDecimal longTermMean;
    private final BigDecimal volatility;
    private final BigDecimal initialRate;
    private final FellerMode fellerMode;

    private CirParameters(BigDecimal meanReversionSpeed, BigDecimal longTermMean, BigDecimal volatility,
                          BigDecimal initialRate, FellerMode fellerMode) {
        this.meanReversionSpeed = requirePositive(meanReversionSpeed, "meanReversionSpeed (theta)");
        this.longTermMean = requirePositive(longTermMean, "longTermMean (mu)");
        this.volatility = requirePositive(volatility, "volatility (sigma)");
        this.initialRate = requirePositive(initialRate, "initialRate (r0)");
        this.fellerMode = Objects.requireNonNull(fellerMode, "fellerMode");
        if (fellerMode == FellerMode.STRICT && !satisfiesFellerCondition()) {
            throw new FellerConditionViolationException(
                    "Feller condition 2*theta*mu >= sigma^2 violated: 2*" + meanReversionSpeed + "*" + longTermMean
                            + " = " + fellerLeftHandSide() + " < " + fellerRightHandSide()
                            + ". Use FellerMode.ALLOW_VIOLATION for explicit robustness experiments.");
        }
    }

    /** Creates parameters in {@link FellerMode#STRICT} mode from decimal strings. */
    public static CirParameters of(String meanReversionSpeed, String longTermMean, String volatility, String initialRate) {
        return of(meanReversionSpeed, longTermMean, volatility, initialRate, FellerMode.STRICT);
    }

    public static CirParameters of(String meanReversionSpeed, String longTermMean, String volatility,
                                   String initialRate, FellerMode mode) {
        return new CirParameters(new BigDecimal(meanReversionSpeed.trim()), new BigDecimal(longTermMean.trim()),
                new BigDecimal(volatility.trim()), new BigDecimal(initialRate.trim()), mode);
    }

    public static CirParameters of(BigDecimal meanReversionSpeed, BigDecimal longTermMean, BigDecimal volatility,
                                   BigDecimal initialRate, FellerMode mode) {
        return new CirParameters(meanReversionSpeed, longTermMean, volatility, initialRate, mode);
    }

    private static BigDecimal requirePositive(BigDecimal value, String name) {
        Objects.requireNonNull(value, name);
        if (value.signum() <= 0) {
            throw new IllegalArgumentException(name + " must be > 0 but was " + value);
        }
        return value;
    }

    /** Exact decimal value of {@code 2 * theta * mu}. */
    public BigDecimal fellerLeftHandSide() {
        return BigDecimal.TWO.multiply(meanReversionSpeed).multiply(longTermMean);
    }

    /** Exact decimal value of {@code sigma^2}. */
    public BigDecimal fellerRightHandSide() {
        return volatility.multiply(volatility);
    }

    public boolean satisfiesFellerCondition() {
        return fellerLeftHandSide().compareTo(fellerRightHandSide()) >= 0;
    }

    public BigDecimal meanReversionSpeed() { return meanReversionSpeed; }
    public BigDecimal longTermMean() { return longTermMean; }
    public BigDecimal volatility() { return volatility; }
    public BigDecimal initialRate() { return initialRate; }
    public FellerMode fellerMode() { return fellerMode; }

    public double meanReversionSpeedDouble() { return toDouble(meanReversionSpeed); }
    public double longTermMeanDouble() { return toDouble(longTermMean); }
    public double volatilityDouble() { return toDouble(volatility); }
    public double initialRateDouble() { return toDouble(initialRate); }

    /** Correctly rounded conversion of a decimal input to binary64 (via the decimal string). */
    static double toDouble(BigDecimal value) {
        return Double.parseDouble(value.toString());
    }

    @Override
    public String toString() {
        return "CirParameters{theta=" + meanReversionSpeed + ", mu=" + longTermMean + ", sigma=" + volatility
                + ", r0=" + initialRate + ", fellerMode=" + fellerMode + ", feller=" + satisfiesFellerCondition() + "}";
    }
}

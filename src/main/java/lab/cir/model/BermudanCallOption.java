package lab.cir.model;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * Bermudan call on a positive cashflow stream: at each exercise date {@code e_j} the holder may buy the
 * remaining cashflows ({@code t_i > e_j}) for the strike {@code K} (a money amount in the unit of the
 * cashflows). With a single exercise date this is the European option on the same underlying.
 */
public final class BermudanCallOption {

    private final CashflowSchedule underlying;
    private final BigDecimal strike;
    private final BigDecimal[] exerciseDates;
    private final double strikeDouble;
    private final double[] exerciseDatesDouble;

    public BermudanCallOption(CashflowSchedule underlying, BigDecimal strike, List<BigDecimal> exerciseDates) {
        this.underlying = Objects.requireNonNull(underlying, "underlying");
        this.strike = Objects.requireNonNull(strike, "strike");
        if (strike.signum() <= 0) {
            throw new IllegalArgumentException("strike must be > 0");
        }
        if (exerciseDates == null || exerciseDates.isEmpty()) {
            throw new IllegalArgumentException("at least one exercise date required");
        }
        this.exerciseDates = exerciseDates.toArray(BigDecimal[]::new);
        for (int j = 0; j < this.exerciseDates.length; j++) {
            if (this.exerciseDates[j].signum() <= 0) {
                throw new IllegalArgumentException("exercise dates must be > 0");
            }
            if (j > 0 && this.exerciseDates[j].compareTo(this.exerciseDates[j - 1]) <= 0) {
                throw new IllegalArgumentException("exercise dates must be strictly increasing");
            }
        }
        if (this.exerciseDates[this.exerciseDates.length - 1].compareTo(underlying.lastTime()) >= 0) {
            throw new IllegalArgumentException("last exercise date must be before the last cashflow");
        }
        this.strikeDouble = CirParameters.toDouble(strike);
        this.exerciseDatesDouble = Arrays.stream(this.exerciseDates).mapToDouble(CirParameters::toDouble).toArray();
    }

    public static BermudanCallOption of(CashflowSchedule underlying, String strike, String... exerciseDates) {
        return new BermudanCallOption(underlying, new BigDecimal(strike.trim()),
                Arrays.stream(exerciseDates).map(s -> new BigDecimal(s.trim())).toList());
    }

    /** Same contract with cashflows and strike scaled by the same factor. */
    public BermudanCallOption scaled(BigDecimal factor) {
        return new BermudanCallOption(underlying.scaled(factor), strike.multiply(factor), Arrays.asList(exerciseDates));
    }

    /** Same underlying and strike, only the given exercise date (European control case). */
    public BermudanCallOption europeanAt(int exerciseIndex) {
        return new BermudanCallOption(underlying, strike, List.of(exerciseDates[exerciseIndex]));
    }

    public BermudanCallOption withStrike(BigDecimal newStrike) {
        return new BermudanCallOption(underlying, newStrike, Arrays.asList(exerciseDates));
    }

    public CashflowSchedule underlying() { return underlying; }
    public BigDecimal strike() { return strike; }
    public double strikeDouble() { return strikeDouble; }
    public int exerciseDateCount() { return exerciseDates.length; }
    public BigDecimal exerciseDate(int j) { return exerciseDates[j]; }
    public double exerciseDateDouble(int j) { return exerciseDatesDouble[j]; }
    public BigDecimal lastExerciseDate() { return exerciseDates[exerciseDates.length - 1]; }

    @Override
    public String toString() {
        return "BermudanCallOption{K=" + strike.toPlainString() + ", exerciseDates="
                + Arrays.stream(exerciseDates).map(BigDecimal::toPlainString).toList() + ", underlying=" + underlying + "}";
    }
}

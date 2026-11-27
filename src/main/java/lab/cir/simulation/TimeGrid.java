package lab.cir.simulation;

import lab.cir.model.BermudanCallOption;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

/**
 * Uniform simulation grid {@code t_k = k / stepsPerYear}, {@code k = 0..N}, ending at the last exercise date.
 * Every exercise date must lie exactly on the grid (checked in exact decimal arithmetic).
 *
 * <p>{@code normalsPerStep = m > 1} enables Brownian coupling: each step consumes {@code m} fine normals and
 * uses {@code dW = sqrt(dt/m) * (Z_1 + ... + Z_m)}.
 */
public final class TimeGrid {

    private final int stepsPerYear;
    private final int normalsPerStep;
    private final int numberOfSteps;
    private final int[] exerciseStepIndex;

    public TimeGrid(BermudanCallOption option, int stepsPerYear, int normalsPerStep) {
        this(exerciseDates(option), stepsPerYear, normalsPerStep);
    }

    /**
     * Grid ending at the last observation date; states are recorded at every observation date
     * (exercise dates of one contract, or the union of exercise dates of a portfolio).
     */
    public TimeGrid(List<BigDecimal> observationDates, int stepsPerYear, int normalsPerStep) {
        if (stepsPerYear <= 0 || normalsPerStep <= 0) {
            throw new IllegalArgumentException("stepsPerYear and normalsPerStep must be positive");
        }
        if (observationDates.isEmpty()) {
            throw new IllegalArgumentException("at least one observation date required");
        }
        this.stepsPerYear = stepsPerYear;
        this.normalsPerStep = normalsPerStep;
        int m = observationDates.size();
        this.exerciseStepIndex = new int[m];
        BigDecimal spy = BigDecimal.valueOf(stepsPerYear);
        for (int j = 0; j < m; j++) {
            BigDecimal k = observationDates.get(j).multiply(spy);
            if (k.signum() <= 0 || k.stripTrailingZeros().scale() > 0) {
                throw new IllegalArgumentException("observation date " + observationDates.get(j)
                        + " is not a positive point of the grid with " + stepsPerYear + " steps per year");
            }
            exerciseStepIndex[j] = k.intValueExact();
            if (j > 0 && exerciseStepIndex[j] <= exerciseStepIndex[j - 1]) {
                throw new IllegalArgumentException("observation dates must be strictly increasing");
            }
        }
        this.numberOfSteps = exerciseStepIndex[m - 1];
    }

    private static List<BigDecimal> exerciseDates(BermudanCallOption option) {
        List<BigDecimal> d = new ArrayList<>();
        for (int j = 0; j < option.exerciseDateCount(); j++) d.add(option.exerciseDate(j));
        return d;
    }

    public int stepsPerYear() { return stepsPerYear; }
    public int normalsPerStep() { return normalsPerStep; }
    public int numberOfSteps() { return numberOfSteps; }
    public int normalsPerPath() { return numberOfSteps * normalsPerStep; }
    public int exerciseStepIndex(int j) { return exerciseStepIndex[j]; }
    public int exerciseDateCount() { return exerciseStepIndex.length; }

    @Override
    public String toString() {
        return "TimeGrid{stepsPerYear=" + stepsPerYear + ", steps=" + numberOfSteps + ", normalsPerStep=" + normalsPerStep + "}";
    }
}

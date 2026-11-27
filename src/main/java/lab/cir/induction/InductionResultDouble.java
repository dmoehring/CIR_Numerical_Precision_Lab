package lab.cir.induction;

import lab.cir.stats.SampleStatsDouble;

/**
 * Result and full diagnostic trail of a backward induction in double. Arrays are {@code [exerciseIndex][path]}.
 * The decisions stored here are the <em>hypothetical</em> decisions of the backward induction; the realised
 * (first) exercise per path is {@link #exerciseIndex}.
 */
public final class InductionResultDouble {

    public final InductionSettings settings;
    public final TrainingStructure structure;
    /** Trained rules; {@code rules[m-1] == null} (no continuation at the last exercise date). */
    public final DateRuleDouble[] rules;
    /** Selected later exercise index per path before the decision at {@code j} ({@code -1}: none). */
    public final int[][] futureExerciseBefore;
    /** Targets {@code Y_j = D(e_j, tau) E_tau} (0 if no later exercise; 0 at the last date). */
    public final double[][] target;
    /** Estimated continuation value per path ({@code NaN} if the path is not in the estimation set). */
    public final double[][] continuation;
    /** Hypothetical exercise decision at {@code e_j}. */
    public final boolean[][] exercise;
    /** Cluster index per path (-1 if not applicable). */
    public final int[][] cluster;
    /** Final exercise index per path after the induction (in-sample policy), -1 if never exercised. */
    public final int[] exerciseIndex;
    /** Discounted option cashflow at t = 0 per path. */
    public final double[] pathValue;
    public SampleStatsDouble price;

    InductionResultDouble(InductionSettings settings, TrainingStructure structure, int dates, int paths) {
        this.settings = settings;
        this.structure = structure;
        this.rules = new DateRuleDouble[dates];
        this.futureExerciseBefore = new int[dates][];
        this.target = new double[dates][];
        this.continuation = new double[dates][];
        this.exercise = new boolean[dates][];
        this.cluster = new int[dates][];
        this.exerciseIndex = new int[paths];
        this.pathValue = new double[paths];
    }

    public int dates() { return rules.length; }
    public int paths() { return pathValue.length; }
}

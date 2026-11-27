package lab.cir.induction;

import lab.cir.stats.SampleStatsDecimal;

import java.math.BigDecimal;

/** BigDecimal counterpart of {@link InductionResultDouble}; {@code continuation} entries are null if not estimated. */
public final class InductionResultDecimal {

    public final InductionSettings settings;
    public final TrainingStructure structure;
    public final DateRuleDecimal[] rules;
    public final int[][] futureExerciseBefore;
    public final BigDecimal[][] target;
    public final BigDecimal[][] continuation;
    public final boolean[][] exercise;
    public final int[][] cluster;
    public final int[] exerciseIndex;
    public final BigDecimal[] pathValue;
    public SampleStatsDecimal price;

    InductionResultDecimal(InductionSettings settings, TrainingStructure structure, int dates, int paths) {
        this.settings = settings;
        this.structure = structure;
        this.rules = new DateRuleDecimal[dates];
        this.futureExerciseBefore = new int[dates][];
        this.target = new BigDecimal[dates][];
        this.continuation = new BigDecimal[dates][];
        this.exercise = new boolean[dates][];
        this.cluster = new int[dates][];
        this.exerciseIndex = new int[paths];
        this.pathValue = new BigDecimal[paths];
    }

    public int dates() { return rules.length; }
    public int paths() { return pathValue.length; }

    /** Whether the clusters were taken from the double states (Mode A). */
    public boolean sharedStructure() { return structure.source() == TrainingStructure.Source.DOUBLE; }
}

package lab.cir.induction;

import lab.cir.stats.SampleStatsDouble;

/** Forward (out-of-sample) application of a fixed trained policy to independent valuation paths (double). */
public final class ValuationResultDouble {

    /** Continuation value from the fixed rule per date and path (0 at the last date, NaN if undefined). */
    public final double[][] continuation;
    /** Hypothetical decision at every date: {@code E > 0 && E > C}. */
    public final boolean[][] decision;
    public final int[][] cluster;
    /** Realised exercise: first date with a positive decision, -1 if none. */
    public final int[] exerciseIndex;
    public final double[] pathValue;
    public SampleStatsDouble price;

    ValuationResultDouble(int dates, int paths) {
        continuation = new double[dates][paths];
        decision = new boolean[dates][paths];
        cluster = new int[dates][paths];
        exerciseIndex = new int[paths];
        pathValue = new double[paths];
    }
}

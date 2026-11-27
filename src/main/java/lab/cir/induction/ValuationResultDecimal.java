package lab.cir.induction;

import lab.cir.stats.SampleStatsDecimal;

import java.math.BigDecimal;

/** BigDecimal counterpart of {@link ValuationResultDouble}; undefined continuation values are null. */
public final class ValuationResultDecimal {

    public final BigDecimal[][] continuation;
    public final boolean[][] decision;
    public final int[][] cluster;
    public final int[] exerciseIndex;
    public final BigDecimal[] pathValue;
    public SampleStatsDecimal price;

    ValuationResultDecimal(int dates, int paths) {
        continuation = new BigDecimal[dates][paths];
        decision = new boolean[dates][paths];
        cluster = new int[dates][paths];
        exerciseIndex = new int[paths];
        pathValue = new BigDecimal[paths];
    }
}

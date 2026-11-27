package lab.cir.induction;

import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;

import java.math.BigDecimal;
import java.util.Comparator;

/**
 * Decision-independent part of the training: for every exercise date the set of eligible paths and,
 * for the cluster method, the quantile partition and its state-space boundaries.
 *
 * <p>Neither depends on exercise decisions taken during the backward induction (ALL: all paths;
 * IN_THE_MONEY: {@code E_j > 0}), so the structure can be computed upfront. This is what makes the two
 * comparison modes explicit:
 * <ul>
 *   <li>Mode A (shared clusters): the BigDecimal induction uses the structure built from the <b>double</b>
 *       states, i.e. identical eligible path indices and identical clusters;</li>
 *   <li>Mode B (separate clusters): each arithmetic uses the structure built from its own states.</li>
 * </ul>
 * The double induction always uses the double structure, so its result is the same in both modes.
 *
 * <p>Boundaries for new (valuation) states: cluster {@code k} covers {@code (u_{k-1}, u_k]} with
 * {@code u_k} = state of the last training path of cluster {@code k}; states below/above the training range
 * belong to the first/last cluster.
 */
public final class TrainingStructure {

    public enum Source { DOUBLE, DECIMAL }

    private final Source source;
    private final int[][] eligible;
    private final QuantilePartition[] partitions;
    private final double[][] upperBoundsDouble;
    private final BigDecimal[][] upperBoundsDecimal;

    private TrainingStructure(Source source, int[][] eligible, QuantilePartition[] partitions,
                              double[][] upperBoundsDouble, BigDecimal[][] upperBoundsDecimal) {
        this.source = source;
        this.eligible = eligible;
        this.partitions = partitions;
        this.upperBoundsDouble = upperBoundsDouble;
        this.upperBoundsDecimal = upperBoundsDecimal;
    }

    public static TrainingStructure fromDouble(PathDataDouble d, InductionSettings s) {
        int m = d.dates();
        int[][] eligible = new int[m][];
        QuantilePartition[] partitions = new QuantilePartition[m];
        double[][] bounds = new double[m][];
        for (int j = 0; j < m; j++) {
            final double[] e = d.exerciseValue[j];
            eligible[j] = selectPaths(d.paths(), s.pathSelection(), p -> e[p] > 0.0);
            if (s.method() == ContinuationMethod.QUANTILE_CLUSTER && j < m - 1) {
                final double[] state = d.rate[j];
                Comparator<Integer> byState = (a, b) -> Double.compare(state[a], state[b]);
                QuantilePartition qp = QuantilePartition.build(eligible[j], d.paths(), byState,
                        s.clusterCount(), s.minPathsPerCluster());
                partitions[j] = qp;
                bounds[j] = new double[Math.max(0, qp.clusterCount() - 1)];
                for (int k = 0; k + 1 < qp.clusterCount(); k++) {
                    bounds[j][k] = state[qp.sortedPath(qp.clusterEnd(k) - 1)];
                }
            }
        }
        return new TrainingStructure(Source.DOUBLE, eligible, partitions, bounds, null);
    }

    public static TrainingStructure fromDecimal(PathDataDecimal d, InductionSettings s) {
        int m = d.dates();
        int[][] eligible = new int[m][];
        QuantilePartition[] partitions = new QuantilePartition[m];
        BigDecimal[][] bounds = new BigDecimal[m][];
        for (int j = 0; j < m; j++) {
            final BigDecimal[] e = d.exerciseValue[j];
            eligible[j] = selectPaths(d.paths(), s.pathSelection(), p -> e[p].signum() > 0);
            if (s.method() == ContinuationMethod.QUANTILE_CLUSTER && j < m - 1) {
                final BigDecimal[] state = d.rate[j];
                Comparator<Integer> byState = (a, b) -> state[a].compareTo(state[b]);
                QuantilePartition qp = QuantilePartition.build(eligible[j], d.paths(), byState,
                        s.clusterCount(), s.minPathsPerCluster());
                partitions[j] = qp;
                bounds[j] = new BigDecimal[Math.max(0, qp.clusterCount() - 1)];
                for (int k = 0; k + 1 < qp.clusterCount(); k++) {
                    bounds[j][k] = state[qp.sortedPath(qp.clusterEnd(k) - 1)];
                }
            }
        }
        return new TrainingStructure(Source.DECIMAL, eligible, partitions, null, bounds);
    }

    private interface PathPredicate {
        boolean test(int p);
    }

    private static int[] selectPaths(int paths, PathSelection selection, PathPredicate itm) {
        if (selection == PathSelection.ALL) {
            int[] all = new int[paths];
            for (int p = 0; p < paths; p++) all[p] = p;
            return all;
        }
        int count = 0;
        for (int p = 0; p < paths; p++) if (itm.test(p)) count++;
        int[] sel = new int[count];
        int i = 0;
        for (int p = 0; p < paths; p++) if (itm.test(p)) sel[i++] = p;
        return sel;
    }

    public Source source() { return source; }
    public int[] eligible(int j) { return eligible[j]; }
    public QuantilePartition partition(int j) { return partitions[j]; }
    /** Upper cluster bounds (K-1 values) in double; only for {@link Source#DOUBLE}. */
    public double[] upperBoundsDouble(int j) { return upperBoundsDouble == null ? null : upperBoundsDouble[j]; }
    /** Upper cluster bounds in BigDecimal; only for {@link Source#DECIMAL}. */
    public BigDecimal[] upperBoundsDecimal(int j) { return upperBoundsDecimal == null ? null : upperBoundsDecimal[j]; }

    public boolean isEligible(int j, int path) {
        return java.util.Arrays.binarySearch(eligible[j], path) >= 0; // eligible[j] is sorted ascending
    }
}

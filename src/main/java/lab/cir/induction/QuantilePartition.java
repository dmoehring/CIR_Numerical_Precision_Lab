package lab.cir.induction;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Partition of the eligible paths of one exercise date into {@code K} (almost) equally sized quantile clusters.
 *
 * <p>Definition: the {@code n} eligible paths are sorted by {@code (state, pathIndex)} ascending; the path
 * index is the deterministic tie-breaker. Cluster {@code k} ({@code 0 <= k < K}) consists of the sorted
 * positions {@code [floor(k n / K), floor((k+1) n / K))}, so cluster sizes differ by at most one.
 * Paths with identical states may therefore end up in two neighbouring clusters (split by path index).
 *
 * <p>Minimum cluster size: {@code K_eff = max(1, min(K_requested, floor(n / minPathsPerCluster)))};
 * if fewer than {@code minPathsPerCluster} paths are eligible, a single cluster is used.
 * This logic works on indices only and is shared by both arithmetics.
 */
public final class QuantilePartition {

    private final int[] sortedPaths;
    private final int[] clusterStart;
    private final int[] clusterOfPath;
    private final int[] rankOfPath;

    private QuantilePartition(int[] sortedPaths, int[] clusterStart, int totalPaths) {
        this.sortedPaths = sortedPaths;
        this.clusterStart = clusterStart;
        this.clusterOfPath = new int[totalPaths];
        this.rankOfPath = new int[totalPaths];
        Arrays.fill(clusterOfPath, -1);
        Arrays.fill(rankOfPath, -1);
        for (int k = 0; k + 1 < clusterStart.length; k++) {
            for (int pos = clusterStart[k]; pos < clusterStart[k + 1]; pos++) {
                clusterOfPath[sortedPaths[pos]] = k;
                rankOfPath[sortedPaths[pos]] = pos;
            }
        }
    }

    public static int effectiveClusterCount(int eligible, int requested, int minPathsPerCluster) {
        if (eligible == 0) {
            return 0;
        }
        return Math.max(1, Math.min(requested, eligible / minPathsPerCluster));
    }

    /**
     * @param eligible    eligible path indices
     * @param totalPaths  total number of paths (size of the lookup tables)
     * @param stateOrder  order of path indices by state value only (ties are broken here by path index)
     */
    public static QuantilePartition build(int[] eligible, int totalPaths, Comparator<Integer> stateOrder,
                                          int requestedClusters, int minPathsPerCluster) {
        Integer[] boxed = new Integer[eligible.length];
        for (int i = 0; i < eligible.length; i++) {
            boxed[i] = eligible[i];
        }
        Arrays.sort(boxed, stateOrder.thenComparingInt(Integer::intValue));
        int[] sorted = new int[boxed.length];
        for (int i = 0; i < boxed.length; i++) {
            sorted[i] = boxed[i];
        }
        int n = sorted.length;
        int k = effectiveClusterCount(n, requestedClusters, minPathsPerCluster);
        int[] start = new int[k + 1];
        for (int c = 0; c <= k; c++) {
            start[c] = (int) ((long) c * n / k);
        }
        return new QuantilePartition(sorted, k == 0 ? new int[]{0} : start, totalPaths);
    }

    public int clusterCount() { return clusterStart.length - 1; }
    public int clusterStart(int k) { return clusterStart[k]; }
    public int clusterEnd(int k) { return clusterStart[k + 1]; }
    public int clusterSize(int k) { return clusterStart[k + 1] - clusterStart[k]; }
    public int sortedPath(int position) { return sortedPaths[position]; }
    public int eligibleCount() { return sortedPaths.length; }
    /** Cluster of a path, or -1 if the path is not eligible. */
    public int clusterOf(int path) { return clusterOfPath[path]; }
    /** Position of a path in the sorted order, or -1. */
    public int rankOf(int path) { return rankOfPath[path]; }

    /** Distance (in sorted positions) of a path to the nearest internal cluster boundary. */
    public int rankDistanceToBoundary(int path) {
        int rank = rankOfPath[path];
        if (rank < 0 || clusterCount() <= 1) {
            return Integer.MAX_VALUE;
        }
        int best = Integer.MAX_VALUE;
        for (int c = 1; c < clusterCount(); c++) {
            // boundary lies between positions start[c]-1 and start[c]
            int d = rank < clusterStart[c] ? clusterStart[c] - 1 - rank : rank - clusterStart[c];
            best = Math.min(best, d);
        }
        return best;
    }
}

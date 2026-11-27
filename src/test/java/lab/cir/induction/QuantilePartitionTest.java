package lab.cir.induction;

import org.junit.jupiter.api.Test;

import java.util.Comparator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class QuantilePartitionTest {

    private static int[] range(int n) {
        int[] a = new int[n];
        for (int i = 0; i < n; i++) a[i] = i;
        return a;
    }

    @Test
    void clustersAreEquallyOrAlmostEquallySized() {
        int n = 103;
        double[] state = new double[n];
        for (int i = 0; i < n; i++) state[i] = Math.sin(i * 1.7); // arbitrary order
        QuantilePartition q = QuantilePartition.build(range(n), n, (a, b) -> Double.compare(state[a], state[b]), 10, 5);
        assertEquals(10, q.clusterCount());
        int total = 0;
        for (int k = 0; k < 10; k++) {
            assertTrue(q.clusterSize(k) == 10 || q.clusterSize(k) == 11, "size " + q.clusterSize(k));
            total += q.clusterSize(k);
        }
        assertEquals(n, total);
        // sorted by state
        for (int pos = 1; pos < n; pos++) {
            assertTrue(state[q.sortedPath(pos - 1)] <= state[q.sortedPath(pos)]);
        }
    }

    @Test
    void tiesAreBrokenDeterministicallyByPathIndex() {
        double[] state = {0.5, 0.1, 0.5, 0.5, 0.1, 0.5};
        Comparator<Integer> cmp = (a, b) -> Double.compare(state[a], state[b]);
        QuantilePartition q = QuantilePartition.build(new int[]{5, 3, 1, 0, 4, 2}, 6, cmp, 2, 1);
        // sorted: (0.1,1) (0.1,4) (0.5,0) (0.5,2) (0.5,3) (0.5,5); clusters [0,3) [3,6)
        int[] expected = {1, 4, 0, 2, 3, 5};
        for (int i = 0; i < 6; i++) assertEquals(expected[i], q.sortedPath(i));
        assertEquals(0, q.clusterOf(0), "path 0 (state 0.5) lands in cluster 0 by index tie-break");
        assertEquals(1, q.clusterOf(2), "path 2 with the same state lands in cluster 1");
        // the out-of-sample boundary is the state of the last path of cluster 0: identical states belong to cluster 0
        double ub = state[q.sortedPath(q.clusterEnd(0) - 1)];
        assertEquals(0, ClusterRuleDouble.assign(new double[]{ub}, 0.5));
    }

    @Test
    void minimumClusterSizeReducesClusterCount() {
        assertEquals(4, QuantilePartition.effectiveClusterCount(100, 10, 25));
        assertEquals(1, QuantilePartition.effectiveClusterCount(10, 10, 25), "small sample: one cluster");
        assertEquals(0, QuantilePartition.effectiveClusterCount(0, 10, 25));
        assertEquals(10, QuantilePartition.effectiveClusterCount(1000, 10, 20));
    }

    @Test
    void assignmentOfNewStatesIncludingOutOfRange() {
        double[] upper = {0.02, 0.03, 0.05}; // 4 clusters: (-inf,.02] (.02,.03] (.03,.05] (.05,inf)
        assertEquals(0, ClusterRuleDouble.assign(upper, -1.0), "below training range -> first cluster");
        assertEquals(0, ClusterRuleDouble.assign(upper, 0.02), "upper bound is inclusive");
        assertEquals(1, ClusterRuleDouble.assign(upper, 0.0200001));
        assertEquals(2, ClusterRuleDouble.assign(upper, 0.05));
        assertEquals(3, ClusterRuleDouble.assign(upper, 0.5), "above training range -> last cluster");
        assertEquals(0, ClusterRuleDouble.assign(new double[0], 0.3), "single cluster");
    }

    @Test
    void rankDistanceToBoundary() {
        double[] state = new double[10];
        for (int i = 0; i < 10; i++) state[i] = i;
        QuantilePartition q = QuantilePartition.build(range(10), 10, (a, b) -> Double.compare(state[a], state[b]), 2, 1);
        assertEquals(0, q.rankDistanceToBoundary(4));
        assertEquals(0, q.rankDistanceToBoundary(5));
        assertEquals(4, q.rankDistanceToBoundary(0));
        assertEquals(4, q.rankDistanceToBoundary(9));
    }
}

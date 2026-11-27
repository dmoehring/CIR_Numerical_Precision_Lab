package lab.cir.induction;

import java.util.Objects;

/** Settings of the backward induction / continuation estimation (identical for both arithmetics). */
public record InductionSettings(ContinuationMethod method,
                                PathSelection pathSelection,
                                int clusterCount,
                                int minPathsPerCluster,
                                RegressionVariable regressionVariable,
                                SummationMode doubleSummation) {

    public InductionSettings {
        Objects.requireNonNull(method);
        Objects.requireNonNull(pathSelection);
        Objects.requireNonNull(regressionVariable);
        Objects.requireNonNull(doubleSummation);
        if (clusterCount < 1 || minPathsPerCluster < 1) {
            throw new IllegalArgumentException("clusterCount and minPathsPerCluster must be >= 1");
        }
    }

    public static InductionSettings clusters(int clusterCount, int minPathsPerCluster) {
        return new InductionSettings(ContinuationMethod.QUANTILE_CLUSTER, PathSelection.ALL, clusterCount,
                minPathsPerCluster, RegressionVariable.SHORT_RATE, SummationMode.NAIVE);
    }

    public InductionSettings withMethod(ContinuationMethod m) {
        return new InductionSettings(m, pathSelection, clusterCount, minPathsPerCluster, regressionVariable, doubleSummation);
    }

    public InductionSettings withPathSelection(PathSelection s) {
        return new InductionSettings(method, s, clusterCount, minPathsPerCluster, regressionVariable, doubleSummation);
    }

    public InductionSettings withSummation(SummationMode s) {
        return new InductionSettings(method, pathSelection, clusterCount, minPathsPerCluster, regressionVariable, s);
    }

    public InductionSettings withClusterCount(int k) {
        return new InductionSettings(method, pathSelection, k, minPathsPerCluster, regressionVariable, doubleSummation);
    }

    public String label() {
        return method == ContinuationMethod.QUANTILE_CLUSTER
                ? "QUANTILE_CLUSTER(K=" + clusterCount + ", min=" + minPathsPerCluster + ", " + pathSelection + ")"
                : "LSM(1,x,x^2; x=" + regressionVariable + ", " + pathSelection + ")";
    }
}

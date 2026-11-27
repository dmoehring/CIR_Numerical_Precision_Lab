package lab.cir.induction;

import java.util.Objects;

/** Settings of the backward induction / continuation estimation (identical for both arithmetics). */
public record InductionSettings(ContinuationMethod method,
                                PathSelection pathSelection,
                                int clusterCount,
                                int minPathsPerCluster,
                                RegressionVariable regressionVariable,
                                SummationMode doubleSummation,
                                int basisDegree) {

    public InductionSettings {
        Objects.requireNonNull(method);
        Objects.requireNonNull(pathSelection);
        Objects.requireNonNull(regressionVariable);
        Objects.requireNonNull(doubleSummation);
        if (clusterCount < 1 || minPathsPerCluster < 1) {
            throw new IllegalArgumentException("clusterCount and minPathsPerCluster must be >= 1");
        }
        if (basisDegree < 0 || basisDegree > 6) {
            throw new IllegalArgumentException("LSM basis degree must be in 0..6");
        }
    }

    /** LSM basis {@code 1, x, x^2} (degree 2). */
    public InductionSettings(ContinuationMethod method, PathSelection pathSelection, int clusterCount,
                             int minPathsPerCluster, RegressionVariable regressionVariable, SummationMode doubleSummation) {
        this(method, pathSelection, clusterCount, minPathsPerCluster, regressionVariable, doubleSummation, 2);
    }

    public static InductionSettings clusters(int clusterCount, int minPathsPerCluster) {
        return new InductionSettings(ContinuationMethod.QUANTILE_CLUSTER, PathSelection.ALL, clusterCount,
                minPathsPerCluster, RegressionVariable.SHORT_RATE, SummationMode.NAIVE);
    }

    public InductionSettings withMethod(ContinuationMethod m) {
        return new InductionSettings(m, pathSelection, clusterCount, minPathsPerCluster, regressionVariable, doubleSummation,
                basisDegree);
    }

    public InductionSettings withPathSelection(PathSelection s) {
        return new InductionSettings(method, s, clusterCount, minPathsPerCluster, regressionVariable, doubleSummation,
                basisDegree);
    }

    public InductionSettings withSummation(SummationMode s) {
        return new InductionSettings(method, pathSelection, clusterCount, minPathsPerCluster, regressionVariable, s,
                basisDegree);
    }

    public InductionSettings withClusterCount(int k) {
        return new InductionSettings(method, pathSelection, k, minPathsPerCluster, regressionVariable, doubleSummation,
                basisDegree);
    }

    public InductionSettings withBasisDegree(int d) {
        return new InductionSettings(method, pathSelection, clusterCount, minPathsPerCluster, regressionVariable,
                doubleSummation, d);
    }

    public String label() {
        return method == ContinuationMethod.QUANTILE_CLUSTER
                ? "QUANTILE_CLUSTER(K=" + clusterCount + ", min=" + minPathsPerCluster + ", " + pathSelection + ")"
                : "LSM(" + basisLabel() + "; x=" + regressionVariable + ", " + pathSelection + ")";
    }

    private String basisLabel() {
        StringBuilder b = new StringBuilder("1");
        for (int k = 1; k <= basisDegree; k++) b.append(k == 1 ? ",x" : ",x^" + k);
        return b.toString();
    }
}

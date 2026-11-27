package lab.cir.induction;

/** Estimator of the conditional continuation value. */
public enum ContinuationMethod {
    /** Primary method: piecewise constant mean of targets within equally sized short-rate quantile clusters. */
    QUANTILE_CLUSTER,
    /** Longstaff-Schwartz least squares regression on {1, x, x^2}, solved by Householder QR. */
    LEAST_SQUARES
}

package lab.cir.induction;

/**
 * Regression rule {@code C(x) = sum_k b_k x^k}, degree {@code coefficients.length - 1} (default 2). Degree 2 is
 * evaluated as {@code b0 + b1 x + b2 (x x)} (the original order), other degrees by Horner's scheme.
 */
public record RegressionRuleDouble(RegressionVariable variable, double[] coefficients, int pathsUsed)
        implements DateRuleDouble {

    public double predict(double x) {
        if (coefficients.length == 3) {
            return coefficients[0] + coefficients[1] * x + coefficients[2] * (x * x);
        }
        double s = coefficients[coefficients.length - 1];
        for (int k = coefficients.length - 2; k >= 0; k--) s = s * x + coefficients[k];
        return s;
    }

    public double predict(double rate, double bondValue) {
        return predict(variable == RegressionVariable.SHORT_RATE ? rate : bondValue);
    }
}

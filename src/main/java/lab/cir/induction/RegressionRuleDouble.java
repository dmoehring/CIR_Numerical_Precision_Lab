package lab.cir.induction;

/** Regression rule {@code C(x) = b0 + b1 x + b2 x^2}. {@code fallbackConstant} if too few paths. */
public record RegressionRuleDouble(RegressionVariable variable, double[] coefficients, int pathsUsed)
        implements DateRuleDouble {

    public double predict(double x) {
        return coefficients[0] + coefficients[1] * x + coefficients[2] * (x * x);
    }

    public double predict(double rate, double bondValue) {
        return predict(variable == RegressionVariable.SHORT_RATE ? rate : bondValue);
    }
}

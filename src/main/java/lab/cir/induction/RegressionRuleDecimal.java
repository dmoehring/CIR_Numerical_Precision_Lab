package lab.cir.induction;

import java.math.BigDecimal;
import java.math.MathContext;

/** BigDecimal regression rule {@code C(x) = b0 + b1 x + b2 x^2} (same evaluation order as double). */
public record RegressionRuleDecimal(RegressionVariable variable, BigDecimal[] coefficients, int pathsUsed)
        implements DateRuleDecimal {

    public BigDecimal predict(BigDecimal x, MathContext mc) {
        if (coefficients.length != 3) {
            BigDecimal s = coefficients[coefficients.length - 1];
            for (int k = coefficients.length - 2; k >= 0; k--) s = s.multiply(x, mc).add(coefficients[k], mc);
            return s;
        }
        return coefficients[0].add(coefficients[1].multiply(x, mc), mc)
                .add(coefficients[2].multiply(x.multiply(x, mc), mc), mc);
    }

    public BigDecimal predict(BigDecimal rate, BigDecimal bondValue, MathContext mc) {
        return predict(variable == RegressionVariable.SHORT_RATE ? rate : bondValue, mc);
    }
}

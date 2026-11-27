package lab.cir.induction;

import java.math.BigDecimal;
import java.math.MathContext;

/** BigDecimal regression rule {@code C(x) = b0 + b1 x + b2 x^2} (same evaluation order as double). */
public record RegressionRuleDecimal(RegressionVariable variable, BigDecimal[] coefficients, int pathsUsed)
        implements DateRuleDecimal {

    public BigDecimal predict(BigDecimal x, MathContext mc) {
        return coefficients[0].add(coefficients[1].multiply(x, mc), mc)
                .add(coefficients[2].multiply(x.multiply(x, mc), mc), mc);
    }

    public BigDecimal predict(BigDecimal rate, BigDecimal bondValue, MathContext mc) {
        return predict(variable == RegressionVariable.SHORT_RATE ? rate : bondValue, mc);
    }
}

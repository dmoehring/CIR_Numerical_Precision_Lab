package lab.cir.analytic;

import lab.cir.model.CirParameters;
import lab.cir.numeric.DecimalMath;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Analytic CIR zero-coupon bond price in {@code BigDecimal}; same formula and the same order of
 * operations as {@link CirZeroBondDouble}. Every operation is rounded to the working MathContext;
 * exp/log/sqrt/expm1 come from {@link DecimalMath} (big-math with guard digits).
 */
public final class CirZeroBondDecimal {

    public record Coefficients(BigDecimal logA, BigDecimal b) {
        public BigDecimal price(BigDecimal rate, DecimalMath math) {
            MathContext mc = math.mc();
            return math.exp(logA.subtract(b.multiply(rate, mc), mc));
        }
    }

    private static final BigDecimal TWO = BigDecimal.TWO;

    private final DecimalMath math;
    private final MathContext mc;
    private final BigDecimal theta;
    private final BigDecimal h;
    private final BigDecimal exponent;
    private final ZeroBondForm form;

    public CirZeroBondDecimal(CirParameters p, DecimalMath math) {
        this(p, math, ZeroBondForm.LOG_DIFFERENCE);
    }

    public CirZeroBondDecimal(CirParameters p, DecimalMath math, ZeroBondForm form) {
        this.form = form;
        this.math = math;
        this.mc = math.mc();
        this.theta = p.meanReversionSpeed();
        BigDecimal mu = p.longTermMean();
        BigDecimal sigma = p.volatility();
        BigDecimal sigma2 = sigma.multiply(sigma, mc);
        this.h = math.sqrt(theta.multiply(theta, mc).add(TWO.multiply(sigma2, mc), mc));
        this.exponent = TWO.multiply(theta, mc).multiply(mu, mc).divide(sigma2, mc);
    }

    public Coefficients coefficients(BigDecimal tau) {
        if (tau.signum() < 0) {
            throw new IllegalArgumentException("tau must be >= 0");
        }
        if (tau.signum() == 0) {
            return new Coefficients(BigDecimal.ZERO, BigDecimal.ZERO);
        }
        BigDecimal thetaPlusH = theta.add(h, mc);
        BigDecimal em1 = math.expm1(h.multiply(tau, mc));
        BigDecimal denom = TWO.multiply(h, mc).add(thetaPlusH.multiply(em1, mc), mc);
        BigDecimal b = TWO.multiply(em1, mc).divide(denom, mc);
        BigDecimal bracket;
        if (form == ZeroBondForm.LOG_DIFFERENCE) {
            bracket = math.log(TWO.multiply(h, mc))
                    .add(thetaPlusH.multiply(tau, mc).divide(TWO, mc), mc)
                    .subtract(math.log(denom), mc);
        } else {
            BigDecimal g = thetaPlusH.multiply(tau, mc).divide(TWO, mc);
            BigDecimal num = TWO.multiply(h, mc).multiply(math.expm1(g), mc).subtract(thetaPlusH.multiply(em1, mc), mc);
            // log1p(x) = log(1 + x); 1 + x is formed at working precision (no double involved)
            bracket = math.log(BigDecimal.ONE.add(num.divide(denom, mc), mc));
        }
        BigDecimal logA = exponent.multiply(bracket, mc);
        return new Coefficients(logA, b);
    }

    public BigDecimal price(BigDecimal tau, BigDecimal rate) {
        return coefficients(tau).price(rate, math);
    }

    public DecimalMath math() { return math; }
}

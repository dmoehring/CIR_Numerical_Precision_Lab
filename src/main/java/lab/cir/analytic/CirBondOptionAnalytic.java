package lab.cir.analytic;

import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;
import lab.cir.model.CirParameters;

/**
 * Double-based analytic reference prices in the CIR model (labelled as such in all reports; this is <b>not</b>
 * a high-precision BigDecimal solution).
 *
 * <p>Zero-bond call (Cox, Ingersoll, Ross 1985; Brigo &amp; Mercurio 2006, eq. (3.26)), in our notation
 * ({@code theta} = speed, {@code mu} = long-term mean), valuation at t = 0, option expiry T, bond maturity S:
 *
 * <pre>
 *   ZBC = P(0,S) chi2(2 r^ (rho + psi + B(T,S)); 4 theta mu / sigma^2, 2 rho^2 r0 e^{hT} / (rho + psi + B(T,S)))
 *       - X P(0,T) chi2(2 r^ (rho + psi);         4 theta mu / sigma^2, 2 rho^2 r0 e^{hT} / (rho + psi))
 *   rho = 2h / (sigma^2 (e^{hT} - 1)),  psi = (theta + h) / sigma^2,  r^ = ln(A(T,S) / X) / B(T,S)
 * </pre>
 *
 * <p>Coupon-bond option via Jamshidian (1989): in the one-factor CIR model every zero bond price is strictly
 * decreasing in {@code r_T}, hence so is the coupon bond. If {@code K < sum c_i A(T,s_i)} (the bond value at
 * {@code r = 0}) there is a unique critical rate {@code r* > 0} with {@code sum c_i P(T,s_i; r*) = K}; then
 * {@code Call = sum c_i ZBC(0, T, s_i, P(T,s_i; r*))}. Otherwise ({@code K >= bond value at r = 0}) the call
 * can never be exercised with profit because {@code r_T >= 0}: its value is 0 (edge case without root).
 */
public final class CirBondOptionAnalytic {

    private final double theta;
    private final double mu;
    private final double sigma;
    private final double r0;
    private final double h;
    private final CirZeroBondDouble zeroBond;

    public CirBondOptionAnalytic(CirParameters p) {
        this.theta = p.meanReversionSpeedDouble();
        this.mu = p.longTermMeanDouble();
        this.sigma = p.volatilityDouble();
        this.r0 = p.initialRateDouble();
        this.zeroBond = new CirZeroBondDouble(p);
        this.h = zeroBond.h();
    }

    public double zeroBond(double maturity) {
        return zeroBond.price(maturity, r0);
    }

    /** European call on a unit zero bond: expiry T, bond maturity S > T, strike X. */
    public double zeroBondCall(double expiry, double bondMaturity, double strike) {
        if (!(bondMaturity > expiry && expiry > 0)) {
            throw new IllegalArgumentException("need 0 < T < S");
        }
        CirZeroBondDouble.Coefficients ts = zeroBond.coefficients(bondMaturity - expiry);
        double logAts = ts.logA();
        double bts = ts.b();
        double pS = zeroBond(bondMaturity);
        double pT = zeroBond(expiry);
        double logStrike = Math.log(strike);
        if (logStrike >= logAts) {
            return 0.0; // r^ <= 0: the bond can never exceed the strike at T
        }
        double rHat = (logAts - logStrike) / bts;
        double sigma2 = sigma * sigma;
        double rho = 2.0 * h / (sigma2 * Math.expm1(h * expiry));
        double psi = (theta + h) / sigma2;
        double dof = 4.0 * theta * mu / sigma2;
        double scaled = 2.0 * rho * rho * r0 * Math.exp(h * expiry);
        double a1 = rho + psi + bts;
        double a2 = rho + psi;
        double c1 = NoncentralChiSquared.cdf(2.0 * rHat * a1, dof, scaled / a1);
        double c2 = NoncentralChiSquared.cdf(2.0 * rHat * a2, dof, scaled / a2);
        return pS * c1 - strike * pT * c2;
    }

    /** European put on a unit zero bond via put-call parity {@code ZBP = ZBC - P(0,S) + X P(0,T)}. */
    public double zeroBondPut(double expiry, double bondMaturity, double strike) {
        return zeroBondCall(expiry, bondMaturity, strike) - zeroBond(bondMaturity) + strike * zeroBond(expiry);
    }

    /** Critical rate r* of Jamshidian's decomposition, or NaN if K >= bond value at r = 0. */
    public double criticalRate(CashflowSchedule cf, double expiry, double strike) {
        int first = firstAfter(cf, expiry);
        double atZero = bondValueAt(cf, first, expiry, 0.0);
        if (strike >= atZero) {
            return Double.NaN;
        }
        double lo = 0.0;
        double hi = Math.max(0.1, 2 * mu);
        while (bondValueAt(cf, first, expiry, hi) > strike) {
            lo = hi;
            hi *= 2.0;
            if (hi > 1e6) throw new ArithmeticException("no bracket for the critical rate");
        }
        for (int it = 0; it < 300 && hi - lo > 0; it++) {
            double mid = 0.5 * (lo + hi);
            if (mid <= lo || mid >= hi) break;
            if (bondValueAt(cf, first, expiry, mid) > strike) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /** European call on the remaining cashflows ({@code t_i > T}) via Jamshidian's decomposition. */
    public double couponBondCall(CashflowSchedule cf, double expiry, double strike) {
        double rStar = criticalRate(cf, expiry, strike);
        if (Double.isNaN(rStar)) {
            return 0.0;
        }
        int first = firstAfter(cf, expiry);
        double sum = 0.0;
        for (int i = first; i < cf.size(); i++) {
            double s = cf.timeDouble(i);
            double ki = zeroBond.price(s - expiry, rStar);
            sum += cf.amountDouble(i) * zeroBondCall(expiry, s, ki);
        }
        return sum;
    }

    /** European reference for a single-exercise-date contract. */
    public double europeanCall(BermudanCallOption option) {
        if (option.exerciseDateCount() != 1) {
            throw new IllegalArgumentException("only single exercise date contracts have an analytic price here");
        }
        return couponBondCall(option.underlying(), option.exerciseDateDouble(0), option.strikeDouble());
    }

    private double bondValueAt(CashflowSchedule cf, int first, double expiry, double rate) {
        double sum = 0.0;
        for (int i = first; i < cf.size(); i++) {
            sum += cf.amountDouble(i) * zeroBond.price(cf.timeDouble(i) - expiry, rate);
        }
        return sum;
    }

    private static int firstAfter(CashflowSchedule cf, double t) {
        for (int i = 0; i < cf.size(); i++) {
            if (cf.timeDouble(i) > t) return i;
        }
        return cf.size();
    }
}

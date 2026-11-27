package lab.cir.numeric;

import java.math.BigDecimal;
import java.math.MathContext;

/** BigDecimal counterpart of {@link HouseholderQrDouble}: same algorithm and order of operations. */
public final class HouseholderQrDecimal {

    private HouseholderQrDecimal() {
    }

    public static BigDecimal[] solve(BigDecimal[][] x, BigDecimal[] y, DecimalMath math) {
        MathContext mc = math.mc();
        int k = x.length;
        int n = y.length;
        if (n < k) {
            throw new IllegalArgumentException("need at least as many rows as columns");
        }
        BigDecimal[] diagonal = new BigDecimal[k];
        for (int c = 0; c < k; c++) {
            BigDecimal[] col = x[c];
            BigDecimal norm2 = BigDecimal.ZERO;
            for (int i = c; i < n; i++) {
                norm2 = norm2.add(col[i].multiply(col[i], mc), mc);
            }
            BigDecimal norm = math.sqrt(norm2);
            if (norm.signum() == 0) {
                diagonal[c] = BigDecimal.ZERO;
                continue;
            }
            BigDecimal alpha = col[c].signum() > 0 ? norm.negate() : norm;
            col[c] = col[c].subtract(alpha, mc);
            BigDecimal vtv = BigDecimal.ZERO;
            for (int i = c; i < n; i++) {
                vtv = vtv.add(col[i].multiply(col[i], mc), mc);
            }
            BigDecimal beta = BigDecimal.TWO.divide(vtv, mc);
            for (int c2 = c + 1; c2 < k; c2++) {
                applyReflector(col, x[c2], c, n, beta, mc);
            }
            applyReflector(col, y, c, n, beta, mc);
            diagonal[c] = alpha;
        }
        BigDecimal[] coefficients = new BigDecimal[k];
        for (int c = k - 1; c >= 0; c--) {
            if (diagonal[c].signum() == 0) {
                coefficients[c] = BigDecimal.ZERO;
                continue;
            }
            BigDecimal s = y[c];
            for (int c2 = c + 1; c2 < k; c2++) {
                s = s.subtract(x[c2][c].multiply(coefficients[c2], mc), mc);
            }
            coefficients[c] = s.divide(diagonal[c], mc);
        }
        return coefficients;
    }

    private static void applyReflector(BigDecimal[] v, BigDecimal[] target, int from, int n, BigDecimal beta, MathContext mc) {
        BigDecimal w = BigDecimal.ZERO;
        for (int i = from; i < n; i++) {
            w = w.add(v[i].multiply(target[i], mc), mc);
        }
        BigDecimal f = beta.multiply(w, mc);
        for (int i = from; i < n; i++) {
            target[i] = target[i].subtract(f.multiply(v[i], mc), mc);
        }
    }
}

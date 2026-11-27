package lab.cir.numeric;

/**
 * Least squares {@code min ||X beta - y||} via Householder QR (no normal equations, no matrix inversion).
 * Column-major input {@code X[col][row]}; the input arrays are overwritten. A column whose remaining norm
 * is zero (rank deficiency) gets coefficient 0.
 */
public final class HouseholderQrDouble {

    private HouseholderQrDouble() {
    }

    public static double[] solve(double[][] x, double[] y) {
        int k = x.length;
        int n = y.length;
        if (n < k) {
            throw new IllegalArgumentException("need at least as many rows as columns");
        }
        double[] diagonal = new double[k];
        for (int c = 0; c < k; c++) {
            double[] col = x[c];
            double norm2 = 0.0;
            for (int i = c; i < n; i++) {
                norm2 += col[i] * col[i];
            }
            double norm = Math.sqrt(norm2);
            if (norm == 0.0) {
                diagonal[c] = 0.0;
                continue;
            }
            double alpha = col[c] > 0.0 ? -norm : norm;
            // v = col[c..n) - alpha e_1, stored in place
            col[c] = col[c] - alpha;
            double vtv = 0.0;
            for (int i = c; i < n; i++) {
                vtv += col[i] * col[i];
            }
            double beta = 2.0 / vtv;
            for (int c2 = c + 1; c2 < k; c2++) {
                applyReflector(col, x[c2], c, n, beta);
            }
            applyReflector(col, y, c, n, beta);
            diagonal[c] = alpha;
        }
        double[] coefficients = new double[k];
        for (int c = k - 1; c >= 0; c--) {
            if (diagonal[c] == 0.0) {
                coefficients[c] = 0.0;
                continue;
            }
            double s = y[c];
            for (int c2 = c + 1; c2 < k; c2++) {
                s -= x[c2][c] * coefficients[c2];
            }
            coefficients[c] = s / diagonal[c];
        }
        return coefficients;
    }

    private static void applyReflector(double[] v, double[] target, int from, int n, double beta) {
        double w = 0.0;
        for (int i = from; i < n; i++) {
            w += v[i] * target[i];
        }
        double f = beta * w;
        for (int i = from; i < n; i++) {
            target[i] -= f * v[i];
        }
    }
}

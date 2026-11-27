package lab.cir.report;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/** Console formatting helpers (output rounding only; never used for computations). */
public final class Fmt {

    private static final MathContext SHOW = new MathContext(18, RoundingMode.HALF_EVEN);

    private Fmt() {
    }

    public static String g(double v) {
        if (Double.isNaN(v)) return "NaN";
        if (v == 0.0) return "0";
        double a = Math.abs(v);
        if (a >= 1e-3 && a < 1e6) return String.format("%.6f", v);
        return String.format("%.3e", v);
    }

    public static String e(double v) {
        if (Double.isNaN(v)) return "NaN";
        if (v == 0.0) return "0";
        return String.format("%.2e", v);
    }

    public static String d(BigDecimal v) {
        if (v == null) return "n/a";
        if (v.signum() == 0) return "0";
        return v.round(SHOW).toString();
    }

    public static String de(BigDecimal v) {
        if (v == null) return "n/a";
        if (v.signum() == 0) return "0";
        return e(v.doubleValue());
    }

    public static String ms(long nanos) {
        return String.format("%.1f", nanos / 1e6);
    }

    public static String mb(long bytes) {
        return bytes < 0 ? "n/a" : String.format("%.1f", bytes / (1024.0 * 1024.0));
    }

    public static String pct(double fraction) {
        return String.format("%.1f%%", 100.0 * fraction);
    }
}

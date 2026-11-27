package lab.cir.metrics;

import lab.cir.numeric.DecimalMath;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.stats.DifferenceStats;
import lab.cir.stats.PairedDifference;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * double vs BigDecimal comparison of two policies evaluated on the same paths (common random numbers):
 * price difference (12.1), paired path-wise analysis (12.11), decision differences (12.6), distance of
 * differing paths to the estimated exercise boundary (12.7), cluster assignment differences (12.5),
 * continuation value differences, approximate exercise boundaries in the short rate (12.9).
 */
public final class PolicyComparison {

    /** Per-date figures. Distances use each variant's own {@code |E - C|}. */
    public record DateFigures(int date, int exerciseDouble, int exerciseDecimal, int realisedDouble, int realisedDecimal,
                              int agree, int differ, int doubleExDecimalCont, int decimalExDoubleCont,
                              Summary distanceDouble, Summary distanceDecimal,
                              int sameCluster, int switchedCluster, DifferenceStats continuationDiff,
                              Boundary boundaryDouble, Boundary boundaryDecimal) {
    }

    /**
     * Approximate exercise boundary in the short rate. For a call on a bond (value decreasing in r) exercise should
     * happen for low rates: {@code maxExercised <= minContinuedItm} if the estimated exercise region is an interval.
     * {@code overlap} counts exercised paths with a rate above the lowest continued in-the-money rate.
     */
    public record Boundary(double maxExercisedRate, double minContinuedItmRate, int overlap) {
        public double midpoint() {
            if (Double.isNaN(maxExercisedRate) || Double.isNaN(minContinuedItmRate)) return Double.NaN;
            return 0.5 * (maxExercisedRate + minContinuedItmRate);
        }
    }

    public final String kind;
    public final BigDecimal priceDouble;
    public final BigDecimal priceDecimal;
    public final BigDecimal priceDifference;
    public final double standardErrorDouble;
    public final BigDecimal standardErrorDecimal;
    public final PairedDifference paired;
    public final DateFigures[] dates;

    public PolicyComparison(PolicyView.OfDouble a, PolicyView.OfDecimal b, PathDataDouble dd, PathDataDecimal dc,
                            DecimalMath math) {
        MathContext mc = math.mc();
        this.kind = a.kind();
        this.priceDouble = new BigDecimal(a.price().mean());
        this.priceDecimal = b.price().mean();
        this.priceDifference = priceDecimal.subtract(priceDouble, mc);
        this.standardErrorDouble = a.price().standardError();
        this.standardErrorDecimal = b.price().standardError();
        this.paired = PairedDifference.of(b.pathValue(), a.pathValue(), b.exerciseIndex(), a.exerciseIndex(), math);

        int m = dd.dates();
        int paths = dd.paths();
        this.dates = new DateFigures[m];
        double[] distD = new double[paths];
        double[] distB = new double[paths];
        for (int j = 0; j < m; j++) {
            int exD = 0, exB = 0, realD = 0, realB = 0, agree = 0, dExBCont = 0, bExDCont = 0, same = 0, switched = 0;
            int nDiff = 0;
            DifferenceStats contDiff = new DifferenceStats(1e-12, mc);
            for (int p = 0; p < paths; p++) {
                boolean xd = a.decision()[j][p];
                boolean xb = b.decision()[j][p];
                if (xd) exD++;
                if (xb) exB++;
                if (a.exerciseIndex()[p] == j) realD++;
                if (b.exerciseIndex()[p] == j) realB++;
                if (xd == xb) {
                    agree++;
                } else {
                    if (xd) dExBCont++; else bExDCont++;
                    double cD = a.continuation()[j][p];
                    BigDecimal cB = b.continuation()[j][p];
                    distD[nDiff] = Double.isNaN(cD) ? Double.NaN : Math.abs(dd.exerciseValue[j][p] - cD);
                    distB[nDiff] = cB == null ? Double.NaN : dc.exerciseValue[j][p].subtract(cB, mc).abs().doubleValue();
                    nDiff++;
                }
                int kd = a.cluster()[j][p];
                int kb = b.cluster()[j][p];
                if (kd >= 0 && kb >= 0) {
                    if (kd == kb) same++; else switched++;
                }
                if (!Double.isNaN(a.continuation()[j][p]) && b.continuation()[j][p] != null) {
                    contDiff.add(b.continuation()[j][p], a.continuation()[j][p]);
                }
            }
            dates[j] = new DateFigures(j, exD, exB, realD, realB, agree, paths - agree, dExBCont, bExDCont,
                    Summary.of(distD, nDiff), Summary.of(distB, nDiff), same, switched, contDiff,
                    boundaryDouble(dd, a.decision()[j], j), boundaryDecimal(dc, b.decision()[j], j));
        }
    }

    static Boundary boundaryDouble(PathDataDouble d, boolean[] decision, int j) {
        double maxEx = Double.NaN;
        double minCont = Double.NaN;
        for (int p = 0; p < d.paths(); p++) {
            double r = d.rate[j][p];
            if (decision[p]) {
                if (Double.isNaN(maxEx) || r > maxEx) maxEx = r;
            } else if (d.exerciseValue[j][p] > 0.0) {
                if (Double.isNaN(minCont) || r < minCont) minCont = r;
            }
        }
        int overlap = 0;
        if (!Double.isNaN(minCont)) {
            for (int p = 0; p < d.paths(); p++) {
                if (decision[p] && d.rate[j][p] > minCont) overlap++;
            }
        }
        return new Boundary(maxEx, minCont, overlap);
    }

    static Boundary boundaryDecimal(PathDataDecimal d, boolean[] decision, int j) {
        BigDecimal maxEx = null;
        BigDecimal minCont = null;
        for (int p = 0; p < d.paths(); p++) {
            BigDecimal r = d.rate[j][p];
            if (decision[p]) {
                if (maxEx == null || r.compareTo(maxEx) > 0) maxEx = r;
            } else if (d.exerciseValue[j][p].signum() > 0) {
                if (minCont == null || r.compareTo(minCont) < 0) minCont = r;
            }
        }
        int overlap = 0;
        if (minCont != null) {
            for (int p = 0; p < d.paths(); p++) {
                if (decision[p] && d.rate[j][p].compareTo(minCont) > 0) overlap++;
            }
        }
        // converted for reporting only (boundary positions are diagnostics, ~1e-16 relative is irrelevant here)
        return new Boundary(maxEx == null ? Double.NaN : maxEx.doubleValue(),
                minCont == null ? Double.NaN : minCont.doubleValue(), overlap);
    }

    public double relativePriceDifference() {
        return priceDouble.signum() == 0 ? Double.NaN : priceDifference.doubleValue() / priceDouble.doubleValue();
    }

    public double differenceOverStandardError() {
        return priceDifference.doubleValue() / standardErrorDouble;
    }

    public int totalDecisionDifferences() {
        int s = 0;
        for (DateFigures f : dates) s += f.differ();
        return s;
    }

    public int totalClusterSwitches() {
        int s = 0;
        for (DateFigures f : dates) s += f.switchedCluster();
        return s;
    }
}

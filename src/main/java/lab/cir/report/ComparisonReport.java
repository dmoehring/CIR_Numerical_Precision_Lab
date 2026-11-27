package lab.cir.report;

import lab.cir.config.ExperimentConfig;
import lab.cir.experiment.ComparisonExperiment;
import lab.cir.experiment.ComparisonExperiment.MethodRun;
import lab.cir.induction.ContinuationMethod;
import lab.cir.metrics.InductionDetails;
import lab.cir.metrics.PolicyComparison;
import lab.cir.metrics.StateDifferences;
import lab.cir.stats.DifferenceStats;
import lab.cir.stats.PairedDifference;
import lab.cir.stats.SampleStatsDecimal;
import lab.cir.stats.SampleStatsDouble;

import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.MathContext;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Console summary and CSV export of a {@link ComparisonExperiment.Run}. */
public final class ComparisonReport {

    private ComparisonReport() {
    }

    public static void print(ComparisonExperiment.Run run, String experimentId, PrintStream out) {
        ExperimentConfig c = run.config;
        MathContext mc = run.setup.math.mc();
        out.println("=====================================================================================");
        out.println(" CIR Numerical Precision Lab - double vs BigDecimal comparison   id=" + experimentId);
        out.println("=====================================================================================");
        out.println("Model     : " + c.model);
        out.println("Feller    : 2*theta*mu = " + c.model.fellerLeftHandSide().toPlainString() + "  sigma^2 = "
                + c.model.fellerRightHandSide().toPlainString());
        out.println("Contract  : " + c.option);
        out.println("Grid      : " + run.setup.grid + "  (Euler full truncation, left rectangle discounting)");
        out.println("Paths     : training " + c.trainingPaths + " (seed " + c.trainingSeed + "), valuation "
                + c.valuationPaths + " (seed " + c.valuationSeed + "), threads " + c.threads);
        out.println("BigDecimal: " + c.precision);
        out.println("Euler     : negative internal states  double " + run.trainingDouble.value().totalNegativeStates()
                + " on " + run.trainingDouble.value().pathsWithNegativeStates() + " paths, BigDecimal "
                + run.trainingDecimal.value().totalNegativeStates() + " on "
                + run.trainingDecimal.value().pathsWithNegativeStates() + " paths (training set)");
        out.println();

        out.println("-- Timing (orientation only: System.nanoTime, single run, no warm-up; allocation = calling thread)");
        out.printf("   %-42s %12s %12s %8s %12s %12s%n", "step", "double ms", "decimal ms", "factor", "dbl alloc MB", "dec alloc MB");
        timing(out, "training paths (simulation + exercise values)", run.trainingDouble, run.trainingDecimal);
        timing(out, "valuation paths (simulation + exercise values)", run.valuationDouble, run.valuationDecimal);
        for (MethodRun m : run.methods) {
            String n = m.method() == ContinuationMethod.QUANTILE_CLUSTER ? "cluster" : "LSM";
            timing(out, n + ": training structure (sort/partition)", m.structureDouble, m.structureDecimal);
            timing(out, n + ": backward induction (dbl vs dec Mode A)", m.inductionDouble, m.inductionDecimalA);
            timing(out, n + ": backward induction (dbl vs dec Mode B)", m.inductionDouble, m.inductionDecimalB);
            timing(out, n + ": out-of-sample valuation (Mode A)", m.valuationDouble, m.valuationDecimalA);
        }
        out.println();

        StateDifferences sd = run.trainingStates;
        out.println("-- State differences at exercise dates, training paths (BigDecimal - double; ULP of the double value)");
        out.printf("   %-5s %10s %10s %9s %10s %10s %10s %10s %10s%n", "date", "max|dr|", "mean|dr|", "maxULP r",
                "max|dI|", "max|dD|", "maxrel D", "max|dB|", "max|dE|");
        for (StateDifferences.DateRow r : sd.dates) {
            out.printf("   %-5s %10s %10s %9.1f %10s %10s %10s %10s %10s%n", c.option.exerciseDate(r.date()).toPlainString(),
                    Fmt.e(r.rate().maxAbs()), Fmt.e(r.rate().meanAbs()), r.rate().maxUlp(), Fmt.e(r.integral().maxAbs()),
                    Fmt.e(r.discount().maxAbs()), Fmt.e(r.discount().maxRel()), Fmt.e(r.bond().maxAbs()),
                    Fmt.e(r.exercise().maxAbs()));
        }
        if (sd.steps.length > 0) {
            out.println("-- Short-rate differences along the grid (" + Math.min(c.tracePaths, c.trainingPaths) + " traced paths)");
            out.printf("   %-6s %10s %10s %10s %10s %10s%n", "t", "max|dr|", "mean|dr|", "maxrel r", "meanULP r", "max|dI|");
            int spy = run.setup.grid.stepsPerYear();
            for (StateDifferences.StepRow r : sd.steps) {
                if (r.step() % Math.max(1, spy / 2) == 0) {
                    out.printf("   %-6.2f %10s %10s %10s %10.2f %10s%n", r.time(), Fmt.e(r.rate().maxAbs()),
                            Fmt.e(r.rate().meanAbs()), Fmt.e(r.rate().maxRel()), r.rate().meanUlp(), Fmt.e(r.integral().maxAbs()));
                }
            }
        }
        out.println();

        for (MethodRun m : run.methods) {
            printMethod(out, run, m, mc);
        }
        if (run.methods.size() > 1) {
            printMethodComparison(out, run, mc);
        }
    }

    private static void timing(PrintStream out, String label, Measured<?> d, Measured<?> b) {
        out.printf("   %-42s %12s %12s %8.1f %12s %12s%n", label, Fmt.ms(d.nanos()), Fmt.ms(b.nanos()),
                d.nanos() == 0 ? Double.NaN : (double) b.nanos() / d.nanos(), Fmt.mb(d.allocatedBytes()), Fmt.mb(b.allocatedBytes()));
    }

    private static void printMethod(PrintStream out, ComparisonExperiment.Run run, MethodRun m, MathContext mc) {
        out.println("== Method: " + m.settings.label() + "   precision " + run.config.precision.precision() + " digits");
        out.printf("   %-26s %22s %11s %22s %11s%n", "variant", "in-sample price", "SE", "out-of-sample price", "SE");
        row(out, "double", m.inductionDouble.value().price, m.valuationDouble.value().price);
        if (m.inductionDoubleKahan != null) {
            row(out, "double + Kahan (control)", m.inductionDoubleKahan.price, m.valuationDoubleKahan.price);
        }
        row(out, "BigDecimal Mode A (shared)", m.inductionDecimalA.value().price, m.valuationDecimalA.value().price);
        row(out, "BigDecimal Mode B (separate)", m.inductionDecimalB.value().price, m.valuationDecimalB.value().price);
        out.println("   price differences V_dec - V_dbl (exact decimal difference):");
        out.printf("   %-24s %12s %11s %10s %12s %12s %24s %10s %8s %8s%n", "comparison", "Delta V", "rel", "|D|/SE",
                "paired mean", "SE_d", "95% CI paired", "max|d_p|", "d_p!=0", "tau diff");
        diff(out, "in-sample  Mode A", m.inSampleA, mc);
        diff(out, "in-sample  Mode B", m.inSampleB, mc);
        diff(out, "out-of-sample Mode A", m.outOfSampleA, mc);
        diff(out, "out-of-sample Mode B", m.outOfSampleB, mc);
        PairedDifference f = m.frozenPolicy;
        out.printf("   %-24s %12s %11s %10s %12s %12s %24s %10s %8d %8s%n", "frozen dbl policy (arith)",
                Fmt.de(f.meanDifference()), "", "", Fmt.de(f.meanDifference()), Fmt.de(f.standardError()),
                "[" + Fmt.de(f.ciLow(mc)) + ", " + Fmt.de(f.ciHigh(mc)) + "]", Fmt.de(f.maxAbsDifference()), f.nonZeroPaths(), "0");

        for (PolicyComparison pc : List.of(m.inSampleA, m.inSampleB, m.outOfSampleA, m.outOfSampleB)) {
            String mode = pc == m.inSampleA || pc == m.outOfSampleA ? "Mode A" : "Mode B";
            out.println("   per exercise date, " + pc.kind + ", " + mode
                    + "  (decisions: hypothetical E>0 && E>C; realised = first exercise)");
            out.printf("     %-5s %8s %8s %8s %8s %7s %7s %7s %11s %11s %11s %8s %10s %21s %21s%n", "date", "ex dbl", "ex dec",
                    "real dbl", "real dec", "differ", "dEx/bC", "bEx/dC", "|E-C| med", "|E-C| mean", "|E-C| max", "cl.sw",
                    "max|dC|", "boundary dbl", "boundary dec");
            for (PolicyComparison.DateFigures d : pc.dates) {
                out.printf("     %-5s %8d %8d %8d %8d %7d %7d %7d %11s %11s %11s %8d %10s %21s %21s%n",
                        run.config.option.exerciseDate(d.date()).toPlainString(), d.exerciseDouble(), d.exerciseDecimal(),
                        d.realisedDouble(), d.realisedDecimal(), d.differ(), d.doubleExDecimalCont(), d.decimalExDoubleCont(),
                        Fmt.e(d.distanceDouble().median()), Fmt.e(d.distanceDouble().mean()), Fmt.e(d.distanceDouble().max()),
                        d.switchedCluster(), Fmt.e(d.continuationDiff().maxAbs()),
                        bracket(d.boundaryDouble()), bracket(d.boundaryDecimal()));
            }
        }
        details(out, run, "Mode A", m.detailsA, m);
        details(out, run, "Mode B", m.detailsB, m);
        out.println();
    }

    private static String bracket(PolicyComparison.Boundary b) {
        return String.format("[%s,%s]%s", Fmt.g(b.maxExercisedRate()), Fmt.g(b.minContinuedItmRate()),
                b.overlap() > 0 ? "x" + b.overlap() : "");
    }

    private static void details(PrintStream out, ComparisonExperiment.Run run, String mode, InductionDetails det, MethodRun m) {
        if (!det.clusters.isEmpty()) {
            double maxAbs = 0, maxOverSe = 0;
            for (InductionDetails.ClusterRow r : det.clusters) {
                if (r.absDifference() != null) maxAbs = Math.max(maxAbs, r.absDifference().doubleValue());
                if (!Double.isNaN(r.differenceOverStandardError())) maxOverSe = Math.max(maxOverSe, r.differenceOverStandardError());
            }
            out.printf("   clusters %s: %d cluster estimates, max |C_dec - C_dbl| = %s, max |dC| / (s/sqrt(n)) = %s%n",
                    mode, det.clusters.size(), Fmt.e(maxAbs), Fmt.e(maxOverSe));
        }
        if (!det.coefficients.isEmpty()) {
            double maxRel = 0;
            for (InductionDetails.CoefficientRow r : det.coefficients) {
                if (!Double.isNaN(r.relDifference())) maxRel = Math.max(maxRel, r.relDifference());
            }
            out.printf("   regression %s: max relative coefficient difference = %s%n", mode, Fmt.e(maxRel));
        }
        if (mode.equals("Mode B") && !det.switches.isEmpty()) {
            StringBuilder sb = new StringBuilder("   cluster switches Mode B (training) per date: ");
            for (InductionDetails.SwitchRow s : det.switches) {
                sb.append(run.config.option.exerciseDate(s.date()).toPlainString()).append(": ").append(s.switched())
                        .append(" (<=1 rank from boundary: ").append(s.switchedWithinOneRank())
                        .append(", max|dr| ").append(Fmt.e(s.maxRateDifference())).append(")  ");
            }
            out.println(sb);
        }
        StringBuilder sb = new StringBuilder("   propagation " + mode + " (date: paths with different later exercise / max|dY|): ");
        for (InductionDetails.PropagationRow p : det.propagation) {
            sb.append(run.config.option.exerciseDate(p.date()).toPlainString()).append(": ").append(p.differentFutureExercise())
                    .append(" / ").append(Fmt.e(p.targetDiff().maxAbs())).append("  ");
        }
        out.println(sb);
    }

    private static void row(PrintStream out, String label, SampleStatsDouble in, SampleStatsDouble oos) {
        out.printf("   %-26s %22s %11s %22s %11s%n", label, Fmt.g(in.mean()), Fmt.e(in.standardError()),
                Fmt.g(oos.mean()), Fmt.e(oos.standardError()));
    }

    private static void row(PrintStream out, String label, SampleStatsDecimal in, SampleStatsDecimal oos) {
        out.printf("   %-26s %22s %11s %22s %11s%n", label, Fmt.d(in.mean()), Fmt.de(in.standardError()),
                Fmt.d(oos.mean()), Fmt.de(oos.standardError()));
    }

    private static void diff(PrintStream out, String label, PolicyComparison pc, MathContext mc) {
        PairedDifference p = pc.paired;
        out.printf("   %-24s %12s %11s %10s %12s %12s %24s %10s %8d %8d%n", label, Fmt.de(pc.priceDifference),
                Fmt.e(pc.relativePriceDifference()), Fmt.e(Math.abs(pc.differenceOverStandardError())),
                Fmt.de(p.meanDifference()), Fmt.de(p.standardError()),
                "[" + Fmt.de(p.ciLow(mc)) + ", " + Fmt.de(p.ciHigh(mc)) + "]", Fmt.de(p.maxAbsDifference()),
                p.nonZeroPaths(), p.pathsWithDifferentExercise());
    }

    private static void printMethodComparison(PrintStream out, ComparisonExperiment.Run run, MathContext mc) {
        MethodRun cl = run.method(ContinuationMethod.QUANTILE_CLUSTER);
        MethodRun ls = run.method(ContinuationMethod.LEAST_SQUARES);
        if (cl == null || ls == null) return;
        out.println("== Method comparison: LSM vs quantile clusters (same simulated paths)");
        BigDecimal inDbl = new BigDecimal(ls.inductionDouble.value().price.mean()).subtract(new BigDecimal(cl.inductionDouble.value().price.mean()), mc);
        BigDecimal outDbl = new BigDecimal(ls.valuationDouble.value().price.mean()).subtract(new BigDecimal(cl.valuationDouble.value().price.mean()), mc);
        BigDecimal inDec = ls.inductionDecimalB.value().price.mean().subtract(cl.inductionDecimalB.value().price.mean(), mc);
        BigDecimal outDec = ls.valuationDecimalB.value().price.mean().subtract(cl.valuationDecimalB.value().price.mean(), mc);
        out.printf("   Delta_method = V_LSM - V_Cluster   double: in-sample %s, out-of-sample %s;  BigDecimal(B): in-sample %s, out-of-sample %s%n",
                Fmt.de(inDbl), Fmt.de(outDbl), Fmt.de(inDec), Fmt.de(outDec));
        out.printf("   Delta_precision (Mode B) = V_dec - V_dbl   cluster: in %s, out %s;   LSM: in %s, out %s%n",
                Fmt.de(cl.inSampleB.priceDifference), Fmt.de(cl.outOfSampleB.priceDifference),
                Fmt.de(ls.inSampleB.priceDifference), Fmt.de(ls.outOfSampleB.priceDifference));
        double ratio = Math.abs(outDbl.doubleValue()) / Math.max(Math.abs(cl.outOfSampleB.priceDifference.doubleValue()), Double.MIN_VALUE);
        out.printf("   |Delta_method| / |Delta_precision| (out-of-sample, double vs Mode B, cluster) = %s%n",
                cl.outOfSampleB.priceDifference.signum() == 0 ? "inf (Delta_precision = 0)" : Fmt.e(ratio));
        double[] rc = ComparisonExperiment.realisedRates(cl.valuationDouble.value().exerciseIndex, run.setup.grid.exerciseDateCount());
        double[] rl = ComparisonExperiment.realisedRates(ls.valuationDouble.value().exerciseIndex, run.setup.grid.exerciseDateCount());
        int paths = cl.valuationDouble.value().exerciseIndex.length;
        int differentRealised = 0;
        for (int p = 0; p < paths; p++) {
            if (cl.valuationDouble.value().exerciseIndex[p] != ls.valuationDouble.value().exerciseIndex[p]) differentRealised++;
        }
        out.printf("   %-6s %14s %14s %22s %22s%n", "date", "real. cluster", "real. LSM", "boundary cluster (dbl)", "boundary LSM (dbl)");
        for (int j = 0; j < rc.length; j++) {
            out.printf("   %-6s %14s %14s %22s %22s%n", run.config.option.exerciseDate(j).toPlainString(), Fmt.pct(rc[j]), Fmt.pct(rl[j]),
                    bracket(cl.outOfSampleA.dates[j].boundaryDouble()), bracket(ls.outOfSampleA.dates[j].boundaryDouble()));
        }
        out.printf("   valuation paths with different realised exercise date (cluster vs LSM, double): %d of %d (%s)%n",
                differentRealised, paths, Fmt.pct((double) differentRealised / paths));
        out.println();
    }

    // ------------------------------------------------------------------------------------------------ CSV

    public static void writeCsv(ComparisonExperiment.Run run, String experimentId, Path dir) throws IOException {
        ExperimentConfig c = run.config;
        MathContext mc = run.setup.math.mc();
        CsvTable summary = new CsvTable("summary", "experiment_id", "method", "variant", "cluster_mode", "sample", "precision",
                "rounding", "paths", "price", "std_error", "ci_low", "ci_high", "train_ms", "valuation_ms");
        CsvTable comparisons = new CsvTable("comparisons", "experiment_id", "method", "cluster_mode", "sample", "price_double",
                "price_decimal", "delta_v", "rel_delta", "abs_delta_over_se", "paired_mean", "paired_se", "paired_ci_low",
                "paired_ci_high", "paired_max_abs", "paired_rmse", "paths_nonzero_diff", "paths_diff_exercise",
                "decision_differences", "cluster_switches");
        CsvTable dates = new CsvTable("exercise_dates", "experiment_id", "method", "cluster_mode", "sample", "date",
                "exercise_double", "exercise_decimal", "realised_double", "realised_decimal", "agree", "differ",
                "double_ex_decimal_cont", "decimal_ex_double_cont", "dist_double_median", "dist_double_mean", "dist_double_max",
                "dist_decimal_median", "dist_decimal_mean", "dist_decimal_max", "same_cluster", "switched_cluster",
                "cont_max_abs_diff", "cont_mean_abs_diff", "cont_max_rel_diff", "boundary_double_max_ex", "boundary_double_min_cont",
                "boundary_double_overlap", "boundary_decimal_max_ex", "boundary_decimal_min_cont", "boundary_decimal_overlap");
        CsvTable clusters = new CsvTable("clusters", "experiment_id", "method", "cluster_mode", "date", "cluster", "size_double",
                "size_decimal", "mean_double", "mean_decimal", "abs_diff", "rel_diff", "min_state", "max_state", "target_sd",
                "se_of_mean", "diff_over_se");
        CsvTable coeffs = new CsvTable("regression_coefficients", "experiment_id", "cluster_mode", "date", "coefficient",
                "double", "decimal", "rel_diff");
        CsvTable prop = new CsvTable("propagation", "experiment_id", "method", "cluster_mode", "date",
                "paths_different_later_exercise", "target_max_abs_diff", "target_mean_abs_diff", "later_decision_differences");
        CsvTable sw = new CsvTable("cluster_switches", "experiment_id", "method", "date", "same", "switched",
                "switched_within_one_rank", "max_rank_distance", "max_rate_diff", "mean_rate_diff");
        CsvTable states = new CsvTable("state_differences", "experiment_id", "date", "quantity", "count", "nonzero",
                "max_abs", "mean_abs", "rmse", "max_rel", "near_zero_refs", "mean_ulp", "max_ulp");
        CsvTable steps = new CsvTable("trajectory_differences", "experiment_id", "step", "time", "rate_max_abs",
                "rate_mean_abs", "rate_max_rel", "rate_mean_ulp", "rate_max_ulp", "integral_max_abs", "integral_mean_abs");
        CsvTable timings = new CsvTable("timings", "experiment_id", "step", "variant", "nanos", "allocated_bytes");

        timeRow(timings, experimentId, "training_paths", "double", run.trainingDouble);
        timeRow(timings, experimentId, "training_paths", "decimal", run.trainingDecimal);
        timeRow(timings, experimentId, "valuation_paths", "double", run.valuationDouble);
        timeRow(timings, experimentId, "valuation_paths", "decimal", run.valuationDecimal);

        for (MethodRun m : run.methods) {
            String meth = m.settings.label();
            String pre = m.method().name();
            timeRow(timings, experimentId, pre + "_structure", "double", m.structureDouble);
            timeRow(timings, experimentId, pre + "_structure", "decimal", m.structureDecimal);
            timeRow(timings, experimentId, pre + "_induction", "double", m.inductionDouble);
            timeRow(timings, experimentId, pre + "_induction_modeA", "decimal", m.inductionDecimalA);
            timeRow(timings, experimentId, pre + "_induction_modeB", "decimal", m.inductionDecimalB);
            timeRow(timings, experimentId, pre + "_valuation", "double", m.valuationDouble);
            timeRow(timings, experimentId, pre + "_valuation_modeA", "decimal", m.valuationDecimalA);
            timeRow(timings, experimentId, pre + "_valuation_modeB", "decimal", m.valuationDecimalB);

            sumRow(summary, experimentId, meth, "double", "-", "in-sample", c, m.inductionDouble.value().price, m.inductionDouble.nanos(), 0);
            sumRow(summary, experimentId, meth, "double", "-", "out-of-sample", c, m.valuationDouble.value().price, m.inductionDouble.nanos(), m.valuationDouble.nanos());
            if (m.inductionDoubleKahan != null) {
                sumRow(summary, experimentId, meth, "double_kahan", "-", "in-sample", c, m.inductionDoubleKahan.price, 0, 0);
                sumRow(summary, experimentId, meth, "double_kahan", "-", "out-of-sample", c, m.valuationDoubleKahan.price, 0, 0);
            }
            sumRow(summary, experimentId, meth, "decimal", "A", "in-sample", c, m.inductionDecimalA.value().price, mc, m.inductionDecimalA.nanos(), 0);
            sumRow(summary, experimentId, meth, "decimal", "B", "in-sample", c, m.inductionDecimalB.value().price, mc, m.inductionDecimalB.nanos(), 0);
            sumRow(summary, experimentId, meth, "decimal", "A", "out-of-sample", c, m.valuationDecimalA.value().price, mc, m.inductionDecimalA.nanos(), m.valuationDecimalA.nanos());
            sumRow(summary, experimentId, meth, "decimal", "B", "out-of-sample", c, m.valuationDecimalB.value().price, mc, m.inductionDecimalB.nanos(), m.valuationDecimalB.nanos());
            sumRow(summary, experimentId, meth, "decimal_frozen_double_policy", "-", "out-of-sample", c, m.frozenPolicyPrice, mc, 0, 0);

            Map<String, PolicyComparison> pcs = new LinkedHashMap<>();
            pcs.put("A|in-sample", m.inSampleA);
            pcs.put("B|in-sample", m.inSampleB);
            pcs.put("A|out-of-sample", m.outOfSampleA);
            pcs.put("B|out-of-sample", m.outOfSampleB);
            for (Map.Entry<String, PolicyComparison> e : pcs.entrySet()) {
                String mode = e.getKey().split("\\|")[0];
                String sample = e.getKey().split("\\|")[1];
                PolicyComparison pc = e.getValue();
                PairedDifference p = pc.paired;
                comparisons.row(experimentId, meth, mode, sample, pc.priceDouble, pc.priceDecimal, pc.priceDifference,
                        pc.relativePriceDifference(), Math.abs(pc.differenceOverStandardError()), p.meanDifference(),
                        p.standardError(), p.ciLow(mc), p.ciHigh(mc), p.maxAbsDifference(), p.rmse(), p.nonZeroPaths(),
                        p.pathsWithDifferentExercise(), pc.totalDecisionDifferences(), pc.totalClusterSwitches());
                for (PolicyComparison.DateFigures d : pc.dates) {
                    dates.row(experimentId, meth, mode, sample, c.option.exerciseDate(d.date()), d.exerciseDouble(),
                            d.exerciseDecimal(), d.realisedDouble(), d.realisedDecimal(), d.agree(), d.differ(),
                            d.doubleExDecimalCont(), d.decimalExDoubleCont(), d.distanceDouble().median(),
                            d.distanceDouble().mean(), d.distanceDouble().max(), d.distanceDecimal().median(),
                            d.distanceDecimal().mean(), d.distanceDecimal().max(), d.sameCluster(), d.switchedCluster(),
                            d.continuationDiff().maxAbs(), d.continuationDiff().meanAbs(), d.continuationDiff().maxRel(),
                            d.boundaryDouble().maxExercisedRate(), d.boundaryDouble().minContinuedItmRate(),
                            d.boundaryDouble().overlap(), d.boundaryDecimal().maxExercisedRate(),
                            d.boundaryDecimal().minContinuedItmRate(), d.boundaryDecimal().overlap());
                }
            }
            for (String mode : List.of("A", "B")) {
                InductionDetails det = mode.equals("A") ? m.detailsA : m.detailsB;
                for (InductionDetails.ClusterRow r : det.clusters) {
                    clusters.row(experimentId, meth, mode, c.option.exerciseDate(r.date()), r.cluster(), r.sizeDouble(),
                            r.sizeDecimal(), r.meanDouble(), r.meanDecimal(), r.absDifference(), r.relDifference(),
                            r.minState(), r.maxState(), r.targetStdDev(), r.standardErrorOfMean(), r.differenceOverStandardError());
                }
                for (InductionDetails.CoefficientRow r : det.coefficients) {
                    coeffs.row(experimentId, mode, c.option.exerciseDate(r.date()), r.coefficient(), r.valueDouble(),
                            r.valueDecimal(), r.relDifference());
                }
                for (InductionDetails.PropagationRow r : det.propagation) {
                    prop.row(experimentId, meth, mode, c.option.exerciseDate(r.date()), r.differentFutureExercise(),
                            r.targetDiff().maxAbs(), r.targetDiff().meanAbs(), r.laterDecisionDifferences());
                }
                if (mode.equals("B")) {
                    for (InductionDetails.SwitchRow r : det.switches) {
                        sw.row(experimentId, meth, c.option.exerciseDate(r.date()), r.same(), r.switched(),
                                r.switchedWithinOneRank(), r.maxRankDistance(), r.maxRateDifference(), r.meanRateDifference());
                    }
                }
            }
        }
        for (StateDifferences.DateRow r : run.trainingStates.dates) {
            BigDecimal date = c.option.exerciseDate(r.date());
            stateRow(states, experimentId, date, "rate", r.rate());
            stateRow(states, experimentId, date, "integral", r.integral());
            stateRow(states, experimentId, date, "discount_factor", r.discount());
            stateRow(states, experimentId, date, "bond_value", r.bond());
            stateRow(states, experimentId, date, "exercise_value", r.exercise());
        }
        for (StateDifferences.StepRow r : run.trainingStates.steps) {
            steps.row(experimentId, r.step(), r.time(), r.rate().maxAbs(), r.rate().meanAbs(), r.rate().maxRel(),
                    r.rate().meanUlp(), r.rate().maxUlp(), r.integral().maxAbs(), r.integral().meanAbs());
        }
        List<CsvTable> all = new ArrayList<>(List.of(summary, comparisons, dates, clusters, coeffs, prop, sw, states, steps, timings));
        for (CsvTable t : all) {
            if (t.size() > 0) t.write(dir);
        }
        Map<String, String> extra = new LinkedHashMap<>();
        extra.put("kind", "comparison (diagnostic run; timings are orientation only)");
        extra.put("negative_states_double", String.valueOf(run.trainingDouble.value().totalNegativeStates()));
        extra.put("negative_states_decimal", String.valueOf(run.trainingDecimal.value().totalNegativeStates()));
        RunMetadata.write(dir, experimentId, c.entries(), extra);
    }

    private static void timeRow(CsvTable t, String id, String step, String variant, Measured<?> m) {
        t.row(id, step, variant, m.nanos(), m.allocatedBytes());
    }

    private static void sumRow(CsvTable t, String id, String method, String variant, String mode, String sample,
                               ExperimentConfig c, SampleStatsDouble s, long trainNs, long valNs) {
        t.row(id, method, variant, mode, sample, "binary64", "-", s.n(), s.mean(), s.standardError(), s.ciLow(), s.ciHigh(),
                trainNs / 1e6, valNs / 1e6);
    }

    private static void sumRow(CsvTable t, String id, String method, String variant, String mode, String sample,
                               ExperimentConfig c, SampleStatsDecimal s, MathContext mc, long trainNs, long valNs) {
        t.row(id, method, variant, mode, sample, c.precision.precision(), c.precision.roundingMode(), s.n(), s.mean(),
                s.standardError(), s.ciLow(mc), s.ciHigh(mc), trainNs / 1e6, valNs / 1e6);
    }

    private static void stateRow(CsvTable t, String id, BigDecimal date, String q, DifferenceStats s) {
        t.row(id, date, q, s.count(), s.nonZeroCount(), s.maxAbs(), s.meanAbs(), s.rmse(), s.maxRel(),
                s.nearZeroReferenceCount(), s.meanUlp(), s.maxUlp());
    }
}

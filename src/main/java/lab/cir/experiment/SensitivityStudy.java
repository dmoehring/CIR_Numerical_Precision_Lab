package lab.cir.experiment;

import lab.cir.config.ExperimentConfig;
import lab.cir.report.CsvTable;
import lab.cir.report.Fmt;
import lab.cir.report.RunMetadata;

import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * One-factor sensitivity series: the configured comparison is repeated with one key changed per run
 * (no Cartesian product). For {@code decimal.precision} series the BigDecimal prices are additionally compared with
 * the run of the highest precision (algorithm reference of the same method).
 * For {@code simulation.stepsPerYear} with {@code simulation.couplingBase} set, {@code normalsPerStep = base / steps}.
 * A key of the form {@code a+b} sets several keys to the same value.
 */
public final class SensitivityStudy {

    private SensitivityStudy() {
    }

    private record Row(String value, String method, ComparisonExperiment.MethodRun m, ComparisonExperiment.Run run) {
    }

    public static void run(ExperimentConfig cfg, String key, List<String> values, PrintStream out) throws IOException {
        String id = RunMetadata.newExperimentId("sensitivity-" + key.replace('+', '_'), cfg.label);
        out.println("== Sensitivity series on " + key + " = " + values + "  id=" + id);
        out.println("   all prices out-of-sample unless stated; Delta = V_dec - V_dbl (exact decimal difference)");
        out.printf("   %-8s %-10s %14s %10s %11s %11s %9s %11s %8s %8s %11s %10s %10s %7s%n", "value", "method", "V double",
                "SE", "Delta A", "Delta B", "|dB|/SE", "Delta B in", "dec.diff", "cl.sw", "frozen", "ms dbl", "ms dec", "factor");
        List<Row> rows = new ArrayList<>();
        String base = cfg.get("simulation.couplingBase");
        for (String v : values) {
            ExperimentConfig c = cfg;
            for (String k : key.split("\\+")) {
                c = c.with(k.trim(), v.trim()); // "a+b" sets several keys to the same value (e.g. training and valuation paths)
            }
            if (key.equals("simulation.stepsPerYear") && base != null) {
                int b = Integer.parseInt(base.trim());
                int spy = Integer.parseInt(v.trim());
                if (b % spy != 0) throw new IllegalArgumentException(spy + " does not divide coupling base " + b);
                c = c.with("simulation.normalsPerStep", String.valueOf(b / spy));
            }
            ComparisonExperiment.Run run = ComparisonExperiment.run(c, false);
            for (ComparisonExperiment.MethodRun m : run.methods) {
                rows.add(new Row(v.trim(), m.method().name(), m, run));
                long msD = run.trainingDouble.nanos() + run.valuationDouble.nanos() + m.structureDouble.nanos()
                        + m.inductionDouble.nanos() + m.valuationDouble.nanos();
                long msB = run.trainingDecimal.nanos() + run.valuationDecimal.nanos() + m.structureDecimal.nanos()
                        + m.inductionDecimalB.nanos() + m.valuationDecimalB.nanos();
                out.printf("   %-8s %-10s %14s %10s %11s %11s %9s %11s %8d %8d %11s %10s %10s %7.1f%n", v.trim(),
                        m.method() == lab.cir.induction.ContinuationMethod.QUANTILE_CLUSTER ? "cluster" : "LSM",
                        Fmt.g(m.valuationDouble.value().price.mean()), Fmt.e(m.valuationDouble.value().price.standardError()),
                        Fmt.de(m.outOfSampleA.priceDifference), Fmt.de(m.outOfSampleB.priceDifference),
                        Fmt.e(Math.abs(m.outOfSampleB.differenceOverStandardError())), Fmt.de(m.inSampleB.priceDifference),
                        m.outOfSampleB.totalDecisionDifferences() + m.inSampleB.totalDecisionDifferences(),
                        m.inSampleB.totalClusterSwitches(), Fmt.de(m.frozenPolicy.meanDifference()),
                        Fmt.ms(msD), Fmt.ms(msB), (double) msB / msD);
            }
        }
        CsvTable t = new CsvTable("sensitivity", "experiment_id", "key", "value", "method", "precision", "price_double_in",
                "price_double_out", "se_double_out", "price_decimal_A_out", "price_decimal_B_out", "price_decimal_B_in",
                "delta_A_out", "delta_B_out", "delta_B_in", "abs_delta_B_out_over_se", "decision_diffs_in_B",
                "decision_diffs_out_B", "cluster_switches_in_B", "frozen_policy_mean_diff", "frozen_policy_max_abs",
                "delta_B_in_vs_highest_precision", "delta_B_out_vs_highest_precision", "ms_double_total", "ms_decimal_total_B",
                "ms_double_induction", "ms_decimal_induction_B", "ms_double_paths", "ms_decimal_paths",
                "alloc_double_paths", "alloc_decimal_paths");
        boolean precisionSeries = key.equals("decimal.precision");
        Map<String, Row> highest = new LinkedHashMap<>();
        if (precisionSeries) {
            for (Row r : rows) {
                Row h = highest.get(r.method);
                if (h == null || r.run.config.precision.precision() > h.run.config.precision.precision()) highest.put(r.method, r);
            }
            out.println("   BigDecimal prices vs the highest precision of the series (same method, Mode B):");
        }
        for (Row r : rows) {
            ComparisonExperiment.MethodRun m = r.m;
            BigDecimal vsIn = null;
            BigDecimal vsOut = null;
            if (precisionSeries) {
                Row h = highest.get(r.method);
                MathContext wide = h.run.setup.math.mc();
                vsIn = m.inductionDecimalB.value().price.mean().subtract(h.m.inductionDecimalB.value().price.mean(), wide);
                vsOut = m.valuationDecimalB.value().price.mean().subtract(h.m.valuationDecimalB.value().price.mean(), wide);
                out.printf("     %-6s %-22s in-sample %s   out-of-sample %s%n", r.value, r.method, Fmt.de(vsIn), Fmt.de(vsOut));
            }
            ComparisonExperiment.Run run = r.run;
            long msD = run.trainingDouble.nanos() + run.valuationDouble.nanos() + m.structureDouble.nanos()
                    + m.inductionDouble.nanos() + m.valuationDouble.nanos();
            long msB = run.trainingDecimal.nanos() + run.valuationDecimal.nanos() + m.structureDecimal.nanos()
                    + m.inductionDecimalB.nanos() + m.valuationDecimalB.nanos();
            t.row(id, key, r.value, r.method, run.config.precision.precision(), m.inductionDouble.value().price.mean(),
                    m.valuationDouble.value().price.mean(), m.valuationDouble.value().price.standardError(),
                    m.valuationDecimalA.value().price.mean(), m.valuationDecimalB.value().price.mean(),
                    m.inductionDecimalB.value().price.mean(), m.outOfSampleA.priceDifference, m.outOfSampleB.priceDifference,
                    m.inSampleB.priceDifference, Math.abs(m.outOfSampleB.differenceOverStandardError()),
                    m.inSampleB.totalDecisionDifferences(), m.outOfSampleB.totalDecisionDifferences(),
                    m.inSampleB.totalClusterSwitches(), m.frozenPolicy.meanDifference(), m.frozenPolicy.maxAbsDifference(),
                    vsIn, vsOut, msD / 1e6, msB / 1e6, m.inductionDouble.nanos() / 1e6, m.inductionDecimalB.nanos() / 1e6,
                    (run.trainingDouble.nanos() + run.valuationDouble.nanos()) / 1e6,
                    (run.trainingDecimal.nanos() + run.valuationDecimal.nanos()) / 1e6,
                    run.trainingDouble.allocatedBytes() + run.valuationDouble.allocatedBytes(),
                    run.trainingDecimal.allocatedBytes() + run.valuationDecimal.allocatedBytes());
        }
        if (cfg.writeCsv) {
            java.nio.file.Path dir = cfg.outputDir.resolve(id);
            t.write(dir);
            Map<String, String> extra = new LinkedHashMap<>();
            extra.put("kind", "sensitivity series");
            extra.put("series.key", key);
            extra.put("series.values", String.join(",", values));
            RunMetadata.write(dir, id, cfg.entries(), extra);
            out.println("CSV written to " + dir.toAbsolutePath());
        }
    }
}

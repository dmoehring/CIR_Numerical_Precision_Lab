package lab.cir.experiment;

import lab.cir.analytic.CirBondOptionAnalytic;
import lab.cir.analytic.CirFiniteDifferenceBermudan;
import lab.cir.config.ExperimentConfig;
import lab.cir.induction.BackwardInductionDecimal;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.TrainingStructure;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;
import lab.cir.random.NormalSource;
import lab.cir.report.CsvTable;
import lab.cir.report.Fmt;
import lab.cir.report.Measured;
import lab.cir.report.RunMetadata;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;

import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * European control experiments (single exercise date = the Bermudan engine without continuation estimation):
 * <ol>
 *   <li>call on a unit zero bond: MC (double and BigDecimal) vs closed form (double-based reference);</li>
 *   <li>call on the coupon bond at every configured exercise date: MC vs Jamshidian vs finite differences;</li>
 *   <li>Euler time-step series with coupled Brownian increments for the coupon-bond case (discretisation bias).</li>
 * </ol>
 * The MC error {@code MC - reference} contains Euler bias, rectangle-rule discounting bias and MC noise.
 */
public final class EuropeanValidation {

    private EuropeanValidation() {
    }

    public static void run(ExperimentConfig cfg, PrintStream out) throws IOException {
        String id = RunMetadata.newExperimentId("european", cfg.label);
        CirBondOptionAnalytic analytic = new CirBondOptionAnalytic(cfg.model);
        CirFiniteDifferenceBermudan fd = new CirFiniteDifferenceBermudan(cfg.model);
        CsvTable table = new CsvTable("european_validation", "experiment_id", "case", "expiry", "strike", "steps_per_year",
                "normals_per_step", "paths", "reference_kind", "reference", "fd_reference", "mc_double", "se_double",
                "z_double", "mc_decimal", "se_decimal", "decimal_minus_double", "ms_double", "ms_decimal");
        out.println("== European control cases  id=" + id);
        out.println("   Reference: double-based closed form (CIR 1985 / Jamshidian 1989) - not a high-precision BigDecimal reference.");
        out.println("   z = (MC_double - reference) / SE_double;  MC error = Euler bias + discounting rule bias + MC noise");
        out.printf("   %-30s %6s %8s %14s %14s %14s %10s %7s %14s %11s %9s%n", "case", "T", "K", "reference", "FD", "MC double",
                "SE", "z", "MC decimal", "dec-dbl", "ms dec");

        String zT = getOr(cfg, "european.zcb.expiry", "1.0");
        String zS = getOr(cfg, "european.zcb.maturity", "5.0");
        String zK = getOr(cfg, "european.zcb.strike", "0.85");
        BermudanCallOption zcbCall = BermudanCallOption.of(CashflowSchedule.single(zS, "1"), zK, zT);
        double zRef = analytic.zeroBondCall(Double.parseDouble(zT), Double.parseDouble(zS), Double.parseDouble(zK));
        validate(cfg, out, table, id, "ZBC unit zero bond S=" + zS, zcbCall, zRef, "closed form", cfg.stepsPerYear, cfg.normalsPerStep,
                fd.price(zcbCall, 2000, 1000).price());

        for (int j = 0; j < cfg.option.exerciseDateCount(); j++) {
            BermudanCallOption eu = cfg.option.europeanAt(j);
            double ref = analytic.europeanCall(eu);
            validate(cfg, out, table, id, "coupon bond call (Jamshidian)", eu, ref, "Jamshidian", cfg.stepsPerYear,
                    cfg.normalsPerStep, fd.price(eu, 2000, 1000).price());
        }

        int base = Integer.parseInt(getOr(cfg, "european.couplingBase", "252"));
        int last = cfg.option.exerciseDateCount() - 1;
        BermudanCallOption eu = cfg.option.europeanAt(last);
        double ref = analytic.europeanCall(eu);
        out.println("   Euler time-step series (coupled increments on a " + base + "/year base grid), European at T="
                + eu.exerciseDate(0).toPlainString() + ", double only:");
        for (String s : getOr(cfg, "european.stepSeries", "4,12,36,84,252").split(",")) {
            int spy = Integer.parseInt(s.trim());
            if (base % spy != 0) throw new IllegalArgumentException("steps " + spy + " do not divide the coupling base " + base);
            validateDoubleOnly(cfg, out, table, id, eu, ref, spy, base / spy);
        }
        if (cfg.writeCsv) {
            java.nio.file.Path dir = cfg.outputDir.resolve(id);
            table.write(dir);
            Map<String, String> extra = new LinkedHashMap<>();
            extra.put("kind", "european validation");
            RunMetadata.write(dir, id, cfg.entries(), extra);
            out.println("CSV written to " + dir.toAbsolutePath());
        }
    }

    private static void validate(ExperimentConfig cfg, PrintStream out, CsvTable table, String id, String label,
                                 BermudanCallOption option, double reference, String refKind, int spy, int nps, double fdRef) {
        Setup s = new Setup(cfg.model, option, spy, nps, cfg.precision, cfg.zeroBondForm);
        NormalSource src = s.source(cfg.valuationSeed);
        InductionSettings is = cfg.inductionSettings(ContinuationMethod.QUANTILE_CLUSTER);
        Measured<InductionResultDouble> d = Measured.of(() -> {
            PathDataDouble pd = s.simulateDouble(src, cfg.valuationPaths, cfg.threads);
            return BackwardInductionDouble.run(pd, TrainingStructure.fromDouble(pd, is), is);
        });
        Measured<InductionResultDecimal> b = Measured.of(() -> {
            PathDataDecimal pc = s.simulateDecimal(src, cfg.valuationPaths, cfg.threads);
            return BackwardInductionDecimal.run(pc, TrainingStructure.fromDecimal(pc, is), is, s.math);
        });
        double mc = d.value().price.mean();
        double se = d.value().price.standardError();
        BigDecimal dec = b.value().price.mean();
        BigDecimal diff = dec.subtract(new BigDecimal(mc), s.math.mc());
        double z = (mc - reference) / se;
        out.printf("   %-30s %6s %8s %14s %14s %14s %10s %7.2f %14s %11s %9s%n", label, option.exerciseDate(0).toPlainString(),
                option.strike().toPlainString(), Fmt.g(reference), Fmt.g(fdRef), Fmt.g(mc), Fmt.e(se), z, Fmt.d(dec),
                Fmt.de(diff), Fmt.ms(b.nanos()));
        table.row(id, label, option.exerciseDate(0), option.strike(), spy, nps, cfg.valuationPaths, refKind, reference, fdRef,
                mc, se, z, dec, b.value().price.standardError(), diff, d.nanos() / 1e6, b.nanos() / 1e6);
    }

    private static void validateDoubleOnly(ExperimentConfig cfg, PrintStream out, CsvTable table, String id,
                                           BermudanCallOption option, double reference, int spy, int nps) {
        Setup s = new Setup(cfg.model, option, spy, nps, cfg.precision, cfg.zeroBondForm);
        NormalSource src = s.source(cfg.valuationSeed);
        InductionSettings is = cfg.inductionSettings(ContinuationMethod.QUANTILE_CLUSTER);
        Measured<InductionResultDouble> d = Measured.of(() -> {
            PathDataDouble pd = s.simulateDouble(src, cfg.valuationPaths, cfg.threads);
            return BackwardInductionDouble.run(pd, TrainingStructure.fromDouble(pd, is), is);
        });
        double mc = d.value().price.mean();
        double se = d.value().price.standardError();
        out.printf("     steps/year %4d (normals/step %3d): MC double %s  SE %s  MC - ref %s  z %.2f%n", spy, nps,
                Fmt.g(mc), Fmt.e(se), Fmt.e(mc - reference), (mc - reference) / se);
        table.row(id, "euler series", option.exerciseDate(0), option.strike(), spy, nps, cfg.valuationPaths, "Jamshidian",
                reference, Double.NaN, mc, se, (mc - reference) / se, null, null, null, d.nanos() / 1e6, null);
    }

    static String getOr(ExperimentConfig cfg, String key, String def) {
        String v = cfg.get(key);
        return v == null || v.isBlank() ? def : v.trim();
    }
}

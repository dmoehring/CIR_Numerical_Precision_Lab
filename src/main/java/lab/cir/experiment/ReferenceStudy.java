package lab.cir.experiment;

import lab.cir.analytic.CirBondOptionAnalytic;
import lab.cir.analytic.CirFiniteDifferenceBermudan;
import lab.cir.config.ExperimentConfig;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PolicyValuation;
import lab.cir.induction.TrainingStructure;
import lab.cir.induction.ValuationResultDouble;
import lab.cir.report.CsvTable;
import lab.cir.report.Fmt;
import lab.cir.report.RunMetadata;
import lab.cir.simulation.PathDataDouble;

import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Error budget against an independent deterministic reference:
 * <ul>
 *   <li>finite-difference Bermudan price with grid refinement (continuous-time model, exact policy on the grid);</li>
 *   <li>double Monte Carlo (cluster and LSM, in-sample and out-of-sample) for an Euler time-step series with coupled
 *       Brownian increments;</li>
 *   <li>the precision effect {@code V_dec - V_dbl} of the configured setup for scale.</li>
 * </ul>
 */
public final class ReferenceStudy {

    private ReferenceStudy() {
    }

    public static void run(ExperimentConfig cfg, PrintStream out) throws IOException {
        String id = RunMetadata.newExperimentId("reference", cfg.label);
        CsvTable table = new CsvTable("reference_study", "experiment_id", "kind", "steps_per_year", "normals_per_step",
                "rate_intervals", "fd_steps_per_year", "method", "sample", "price", "std_error", "minus_fd");
        CirFiniteDifferenceBermudan fd = new CirFiniteDifferenceBermudan(cfg.model);
        CirBondOptionAnalytic analytic = new CirBondOptionAnalytic(cfg.model);
        out.println("== Error budget vs deterministic finite-difference reference  id=" + id);
        out.println("   FD grid refinement (Bermudan; European at first date vs Jamshidian as FD sanity check):");
        double fdBest = Double.NaN;
        double previous = Double.NaN;
        int[][] grids = {{250, 100}, {500, 200}, {1000, 400}, {2000, 800}, {4000, 1600}};
        double euRef = analytic.europeanCall(cfg.option.europeanAt(0));
        for (int[] g : grids) {
            CirFiniteDifferenceBermudan.Result r = fd.price(cfg.option, g[0], g[1]);
            double eu = fd.price(cfg.option.europeanAt(0), g[0], g[1]).price();
            out.printf("     rate intervals %5d, steps/year %5d, rMax %.3f:  Bermudan %s  change %s   European(e1) FD-Jamshidian %s%n",
                    g[0], g[1], r.rMax(), Fmt.g(r.price()), Double.isNaN(previous) ? "-" : Fmt.e(r.price() - previous),
                    Fmt.e(eu - euRef));
            table.row(id, "fd", "", "", g[0], g[1], "FD", "-", r.price(), "", "");
            previous = r.price();
            fdBest = r.price();
        }
        out.println("   European at e1 (Jamshidian) = " + Fmt.g(euRef) + "; Bermudan must be >= every European at an exercise date.");

        int base = Integer.parseInt(EuropeanValidation.getOr(cfg, "european.couplingBase", "252"));
        out.println("   Monte Carlo (double) vs FD = " + Fmt.g(fdBest) + "  [training " + cfg.trainingPaths + " / valuation "
                + cfg.valuationPaths + " paths, coupled increments, base " + base + "/year]");
        out.printf("     %-6s %-18s %16s %10s %12s %16s %10s %12s%n", "steps", "method", "in-sample", "SE", "minus FD",
                "out-of-sample", "SE", "minus FD");
        for (String s : EuropeanValidation.getOr(cfg, "european.stepSeries", "4,12,36,84,252").split(",")) {
            int spy = Integer.parseInt(s.trim());
            if (base % spy != 0) throw new IllegalArgumentException("steps " + spy + " do not divide " + base);
            Setup st = new Setup(cfg.model, cfg.option, spy, base / spy, cfg.precision, cfg.zeroBondForm);
            PathDataDouble train = st.simulateDouble(st.source(cfg.trainingSeed), cfg.trainingPaths, cfg.threads);
            PathDataDouble val = st.simulateDouble(st.source(cfg.valuationSeed), cfg.valuationPaths, cfg.threads);
            for (ContinuationMethod m : ContinuationMethod.values()) {
                InductionSettings is = cfg.inductionSettings(m);
                InductionResultDouble ind = BackwardInductionDouble.run(train, TrainingStructure.fromDouble(train, is), is);
                ValuationResultDouble v = PolicyValuation.valueDouble(val, ind);
                out.printf("     %-6d %-18s %16s %10s %12s %16s %10s %12s%n", spy, m == ContinuationMethod.QUANTILE_CLUSTER
                                ? "cluster K=" + cfg.clusterCount : "LSM " + cfg.regressionVariable, Fmt.g(ind.price.mean()),
                        Fmt.e(ind.price.standardError()), Fmt.e(ind.price.mean() - fdBest), Fmt.g(v.price.mean()),
                        Fmt.e(v.price.standardError()), Fmt.e(v.price.mean() - fdBest));
                table.row(id, "mc", spy, base / spy, "", "", m, "in-sample", ind.price.mean(), ind.price.standardError(),
                        ind.price.mean() - fdBest);
                table.row(id, "mc", spy, base / spy, "", "", m, "out-of-sample", v.price.mean(), v.price.standardError(),
                        v.price.mean() - fdBest);
            }
        }
        ComparisonExperiment.Run run = ComparisonExperiment.run(cfg, false);
        out.println("   Precision effect at the configured grid (" + cfg.stepsPerYear + "/year, " + cfg.precision.precision()
                + " digits), out-of-sample, Mode B:");
        for (ComparisonExperiment.MethodRun m : run.methods) {
            out.printf("     %-50s V_dec - V_dbl = %s   (SE %s)%n", m.settings.label(), Fmt.de(m.outOfSampleB.priceDifference),
                    Fmt.e(m.outOfSampleB.standardErrorDouble));
            table.row(id, "precision", cfg.stepsPerYear, cfg.normalsPerStep, "", "", m.method(), "out-of-sample-delta",
                    m.outOfSampleB.priceDifference, m.outOfSampleB.standardErrorDouble, "");
        }
        if (cfg.writeCsv) {
            java.nio.file.Path dir = cfg.outputDir.resolve(id);
            table.write(dir);
            Map<String, String> extra = new LinkedHashMap<>();
            extra.put("kind", "reference study (FD vs MC)");
            RunMetadata.write(dir, id, cfg.entries(), extra);
            out.println("CSV written to " + dir.toAbsolutePath());
        }
    }
}

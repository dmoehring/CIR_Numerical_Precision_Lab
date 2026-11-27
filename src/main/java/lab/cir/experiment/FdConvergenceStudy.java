package lab.cir.experiment;

import lab.cir.analytic.CirBondOptionAnalytic;
import lab.cir.analytic.CirFiniteDifferenceBermudan;
import lab.cir.analytic.CirFiniteDifferenceBermudan.Scheme;
import lab.cir.analytic.CirQuadratureBermudan;
import lab.cir.config.ExperimentConfig;
import lab.cir.report.CsvTable;
import lab.cir.report.RunMetadata;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

/**
 * Accuracy of the two deterministic Bermudan references against each other.
 * <ol>
 *   <li>Quadrature reference ({@link CirQuadratureBermudan}) under refinement of its settings, and its European
 *       prices against Jamshidian's closed form.</li>
 *   <li>Finite differences ({@link CirFiniteDifferenceBermudan}) for four schemes (uniform / sinh grid, with and
 *       without kink smoothing) under joint refinement of rate and time grid: error against the quadrature price,
 *       error ratio per doubling (4 = clean second order) and Richardson extrapolation.</li>
 *   <li>Separate refinement of the rate grid (fine time grid), of the time grid (fine rate grid) and of the
 *       truncation {@code rMax}, for the default scheme.</li>
 *   <li>Exercise boundaries {@code r*_j} of both references.</li>
 * </ol>
 */
public final class FdConvergenceStudy {

    private FdConvergenceStudy() {
    }

    public static void run(ExperimentConfig cfg, PrintStream out) throws IOException {
        String id = RunMetadata.newExperimentId("fdstudy", cfg.label);
        CsvTable table = new CsvTable("fd_study", "experiment_id", "part", "scheme", "rate_intervals", "steps_per_year",
                "rate_max", "price", "minus_quadrature", "error_ratio", "richardson_minus_quadrature");
        out.println("== Deterministic Bermudan references: quadrature vs finite differences  id=" + id);

        // 1. quadrature
        CirQuadratureBermudan quadrature = new CirQuadratureBermudan(cfg.model);
        CirQuadratureBermudan.Settings[] qs = {
                new CirQuadratureBermudan.Settings(1.0, 32, 12, 1e-11),
                CirQuadratureBermudan.Settings.standard(),
                new CirQuadratureBermudan.Settings(1.0, 128, 20, 1e-13),
                new CirQuadratureBermudan.Settings(1.5, 96, 16, 1e-13),
        };
        out.println("   1. Quadrature (exact forward-measure transition, piecewise Chebyshev continuation):");
        CirQuadratureBermudan.Result reference = null;
        for (CirQuadratureBermudan.Settings s : qs) {
            long t0 = System.nanoTime();
            CirQuadratureBermudan.Result r = quadrature.price(cfg.option, s);
            out.printf(Locale.ROOT, "     rMax %.1f panels %4d degree %3d tol %.0e:  %.15f   (%,d density evaluations, %.1f s)%n",
                    s.rateMax(), s.panels(), s.degree(), s.tolerance(), r.price(), r.densityEvaluations(),
                    (System.nanoTime() - t0) / 1e9);
            table.row(id, "quadrature", s.toString().replace(',', ';'), "", "", s.rateMax(), r.price(), "", "", "");
            if (s.equals(CirQuadratureBermudan.Settings.standard())) reference = r;
        }
        double ref = reference.price();
        CirBondOptionAnalytic analytic = new CirBondOptionAnalytic(cfg.model);
        for (int j = 0; j < cfg.option.exerciseDateCount(); j++) {
            double q = quadrature.price(cfg.option.europeanAt(j)).price();
            double jam = analytic.europeanCall(cfg.option.europeanAt(j));
            out.printf(Locale.ROOT, "     European e%d: quadrature %.15f  Jamshidian %.15f  diff %+.2e%n", j + 1, q, jam, q - jam);
            table.row(id, "european", "e" + (j + 1), "", "", "", q, q - jam, "", "");
        }

        // 2. FD schemes, joint refinement
        CirFiniteDifferenceBermudan fd = new CirFiniteDifferenceBermudan(cfg.model);
        Scheme[] schemes = {fd.legacyScheme(), fd.legacyScheme().withSmoothing(true),
                fd.defaultScheme().withSmoothing(false), fd.defaultScheme()};
        int levels = Integer.parseInt(EuropeanValidation.getOr(cfg, "fd.levels", "6"));
        out.println("   2. FD joint refinement (rate intervals x steps/year), error vs quadrature " + String.format(Locale.ROOT, "%.15f", ref));
        for (Scheme s : schemes) {
            out.println("     scheme " + s.label());
            double prev = Double.NaN;
            double prevErr = Double.NaN;
            for (int k = 0; k < levels; k++) {
                int n = 250 << k;
                int spy = 100 << k;
                double p = fd.price(cfg.option, n, spy, s).price();
                double err = p - ref;
                double ratio = prevErr / err;
                double rich = Double.isNaN(prev) ? Double.NaN : p + (p - prev) / 3.0 - ref;
                out.printf(Locale.ROOT, "       %6d x %5d  %.12f  error %+.3e  ratio %6s  Richardson error %s%n", n, spy, p, err,
                        Double.isNaN(ratio) ? "-" : String.format(Locale.ROOT, "%.2f", ratio),
                        Double.isNaN(rich) ? "-" : String.format(Locale.ROOT, "%+.3e", rich));
                table.row(id, "joint", s.label(), n, spy, s.rateMax(), p, err, Double.isNaN(ratio) ? "" : ratio,
                        Double.isNaN(rich) ? "" : rich);
                prev = p;
                prevErr = err;
            }
        }

        // 3. separate refinements, default scheme
        Scheme def = fd.defaultScheme();
        int fineSteps = 6400;
        int fineRate = 16000;
        out.println("   3a. rate grid only (" + fineSteps + " steps/year), " + def.label());
        separate(out, table, id, "space", fd, cfg, def, new int[]{250, 500, 1000, 2000, 4000, 8000}, null, fineSteps, ref);
        out.println("   3b. time grid only (" + fineRate + " rate intervals), " + def.label());
        separate(out, table, id, "time", fd, cfg, def, null, new int[]{25, 50, 100, 200, 400, 800, 1600}, fineRate, ref);
        out.println("   3c. truncation rMax (4000 x 1600), " + def.label());
        for (double rMax : new double[]{0.2, 0.3, 0.5, 1.0}) {
            Scheme s = def.withRateMax(rMax);
            CirFiniteDifferenceBermudan.Result r = fd.price(cfg.option, 4000, 1600, s);
            out.printf(Locale.ROOT, "       rMax %.3f (grid end %.4f): %.12f  error %+.3e%n", rMax, r.rMax(), r.price(), r.price() - ref);
            table.row(id, "rmax", s.label(), 4000, 1600, r.rMax(), r.price(), r.price() - ref, "", "");
        }

        // 4. exercise boundaries
        CirFiniteDifferenceBermudan.Result fine = fd.price(cfg.option, 8000, 3200, def);
        out.println("   4. exercise boundaries r*_j (exercise below; last date: bond value = strike)");
        for (int j = 0; j < cfg.option.exerciseDateCount(); j++) {
            double[] q = reference.exerciseBoundaries()[j];
            double[] f = fine.exerciseBoundaries()[j];
            out.printf(Locale.ROOT, "       e%d: quadrature %s  FD(8000x3200) %s%n", j + 1, Arrays.toString(q), Arrays.toString(f));
            table.row(id, "boundary", "e" + (j + 1), 8000, 3200, "", q.length > 0 ? q[0] : "",
                    q.length > 0 && f.length > 0 ? f[0] - q[0] : "", "", "");
        }

        if (cfg.writeCsv) {
            java.nio.file.Path dir = cfg.outputDir.resolve(id);
            table.write(dir);
            Map<String, String> extra = new LinkedHashMap<>();
            extra.put("kind", "deterministic reference study (quadrature vs FD)");
            RunMetadata.write(dir, id, cfg.entries(), extra);
            out.println("CSV written to " + dir.toAbsolutePath());
        }
    }

    private static void separate(PrintStream out, CsvTable table, String id, String part, CirFiniteDifferenceBermudan fd,
                                 ExperimentConfig cfg, Scheme s, int[] rates, int[] steps, int fixed, double ref) {
        int count = rates != null ? rates.length : steps.length;
        double prevErr = Double.NaN;
        for (int k = 0; k < count; k++) {
            int n = rates != null ? rates[k] : fixed;
            int spy = rates != null ? fixed : steps[k];
            double p = fd.price(cfg.option, n, spy, s).price();
            double err = p - ref;
            double ratio = prevErr / err;
            out.printf(Locale.ROOT, "       %6d x %5d  %.12f  error %+.3e  ratio %s%n", n, spy, p, err,
                    Double.isNaN(ratio) ? "-" : String.format(Locale.ROOT, "%.2f", ratio));
            table.row(id, part, s.label(), n, spy, s.rateMax(), p, err, Double.isNaN(ratio) ? "" : ratio, "");
            prevErr = err;
        }
    }
}

package lab.cir.experiment;

import lab.cir.analytic.CirQuadratureBermudan;
import lab.cir.analytic.CirQuadratureBermudan.ExercisePolicy;
import lab.cir.config.ExperimentConfig;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ClusterRuleDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.DateRuleDouble;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PathSelection;
import lab.cir.induction.RegressionRuleDouble;
import lab.cir.induction.TrainingStructure;
import lab.cir.report.CsvTable;
import lab.cir.report.RunMetadata;
import lab.cir.simulation.ExactPathGenerator;
import lab.cir.simulation.PathDataDouble;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Clean comparison of continuation estimators (quantile clusters vs least squares) against the exact solution.
 *
 * <p>For every training-path source (exact simulation at the exercise dates or Euler paths), path count {@code N}
 * and replication (independent training seed), all estimator configurations are trained on the <b>same</b> paths
 * (paired design). Each trained exercise rule is then valued <b>exactly</b> in the continuous-time model with
 * {@link CirQuadratureBermudan#valuePolicy} - no valuation noise, no Euler bias in the valuation. Metrics:
 * <ul>
 *   <li>policy loss {@code L = V* - V(rule) >= 0} (the only quantity that matters for the price of the policy);</li>
 *   <li>exercise boundary error {@code r^_j - r*_j} (switch points of the rule vs the exact boundary), and the
 *       number of switch points (1 for a threshold rule);</li>
 *   <li>continuation error at the exact boundary {@code C^_j(r*_j) - C*_j(r*_j)} (the decision-relevant error) and
 *       its RMS over the training states;</li>
 *   <li>in-sample price minus {@code V*} (foresight bias plus, for Euler paths, discretisation bias).</li>
 * </ul>
 * Summaries report mean, standard error, median and 90% quantile over replications, and paired differences of every
 * configuration against a baseline configuration.
 */
public final class MethodStudy {

    enum PathSource { EXACT, EULER }

    record Config(String name, InductionSettings settings) {
    }

    record Replication(PathSource source, int paths, int rep, double[] loss, double[][] boundaryError,
                       int[][] switchCount, double[][] continuationErrorAtBoundary, double[] continuationRms,
                       double[] inSampleBias) {
    }

    private MethodStudy() {
    }

    public static void run(ExperimentConfig cfg, PrintStream out) throws Exception {
        String id = RunMetadata.newExperimentId("methodstudy", cfg.label);
        int[] pathCounts = ints(EuropeanValidation.getOr(cfg, "methodstudy.paths", "250,500,1000,2000,4000,8000,16000"));
        int reps = Integer.parseInt(EuropeanValidation.getOr(cfg, "methodstudy.replications", "40"));
        List<PathSource> sources = new ArrayList<>();
        for (String s : EuropeanValidation.getOr(cfg, "methodstudy.sources", "EXACT,EULER").split(",")) {
            sources.add(PathSource.valueOf(s.trim()));
        }
        int threads = Integer.parseInt(EuropeanValidation.getOr(cfg, "methodstudy.threads",
                String.valueOf(Runtime.getRuntime().availableProcessors())));
        List<Config> configs = configs(cfg);
        int baseline = baselineIndex(configs, EuropeanValidation.getOr(cfg, "methodstudy.baseline", "LSM2-ALL"));
        CirQuadratureBermudan.Settings qs = new CirQuadratureBermudan.Settings(1.0,
                Integer.parseInt(EuropeanValidation.getOr(cfg, "methodstudy.quadPanels", "32")),
                Integer.parseInt(EuropeanValidation.getOr(cfg, "methodstudy.quadDegree", "12")), 1e-11);

        CirQuadratureBermudan quad = new CirQuadratureBermudan(cfg.model);
        CirQuadratureBermudan.Result optimal = quad.price(cfg.option, qs);
        double vStar = optimal.price();
        int m = cfg.option.exerciseDateCount();
        int dates = m - 1; // dates with a continuation decision
        double[] rStar = new double[dates];
        for (int j = 0; j < dates; j++) rStar[j] = optimal.exerciseBoundaries()[j][0];

        out.println("== Cluster vs LSM against the exact solution  id=" + id);
        out.printf(Locale.ROOT, "   V* = %.12f (quadrature %d x %d), exact boundaries r* = %s%n", vStar, qs.panels(), qs.degree(),
                Arrays.toString(Arrays.stream(rStar).mapToObj(x -> String.format(Locale.ROOT, "%.6f", x)).toArray()));
        out.println("   sources " + sources + ", paths " + Arrays.toString(pathCounts) + ", " + reps
                + " replications (independent training seeds), all configurations on the same paths; Euler: "
                + cfg.stepsPerYear + " steps/year");
        for (int c = 0; c < configs.size(); c++) {
            out.println("   [" + c + "] " + configs.get(c).name() + " = " + configs.get(c).settings().label()
                    + (c == baseline ? "   (baseline for paired differences)" : ""));
        }

        ExactPathGenerator exact = new ExactPathGenerator(cfg.model, cfg.option);
        Setup euler = Setup.of(cfg);
        long t0 = System.nanoTime();
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        List<Future<Replication>> futures = new ArrayList<>();
        for (PathSource src : sources) {
            for (int n : pathCounts) {
                for (int rep = 0; rep < reps; rep++) {
                    final int r = rep;
                    long seed = cfg.trainingSeed + 7919L * rep;
                    futures.add(pool.submit(() -> replicate(src, n, r, seed, exact, euler, configs, quad, cfg, qs, vStar,
                            optimal, rStar)));
                }
            }
        }
        List<Replication> results = new ArrayList<>();
        for (Future<Replication> f : futures) results.add(f.get());
        pool.shutdown();
        out.printf(Locale.ROOT, "   %d trained policies valued exactly in %.0f s (%d threads)%n",
                results.size() * configs.size(), (System.nanoTime() - t0) / 1e9, threads);

        CsvTable raw = new CsvTable("methodstudy_raw", "experiment_id", "source", "paths", "replication", "config",
                "policy_loss", "in_sample_minus_vstar", "continuation_rms", "boundary_error_e1", "boundary_error_e2",
                "boundary_error_e3", "switch_points_e1", "switch_points_e2", "switch_points_e3",
                "cont_error_at_boundary_e1", "cont_error_at_boundary_e2", "cont_error_at_boundary_e3");
        for (Replication r : results) {
            for (int c = 0; c < configs.size(); c++) {
                List<Object> row = new ArrayList<>(List.of(id, r.source(), r.paths(), r.rep(), configs.get(c).name(),
                        r.loss()[c], r.inSampleBias()[c], r.continuationRms()[c]));
                for (int j = 0; j < 3; j++) row.add(j < dates ? r.boundaryError()[c][j] : "");
                for (int j = 0; j < 3; j++) row.add(j < dates ? r.switchCount()[c][j] : "");
                for (int j = 0; j < 3; j++) row.add(j < dates ? r.continuationErrorAtBoundary()[c][j] : "");
                raw.row(row.toArray());
            }
        }
        CsvTable summary = new CsvTable("methodstudy_summary", "experiment_id", "source", "paths", "config",
                "loss_mean", "loss_se", "loss_median", "loss_p90", "paired_diff_vs_baseline", "paired_diff_se",
                "share_better_than_baseline", "boundary_abs_error_mean", "boundary_error_mean", "non_threshold_share",
                "cont_error_at_boundary_rms", "continuation_rms_mean", "in_sample_bias_mean", "in_sample_bias_se");

        for (PathSource src : sources) {
            out.println();
            out.println("   --- training paths: " + src + " ---  policy loss V* - V(rule) [mean (SE) | median | p90],"
                    + " paired diff vs [" + baseline + "] (SE, share better), mean |r^ - r*|, RMS C-error at r*,"
                    + " in-sample - V*");
            for (int n : pathCounts) {
                List<Replication> rs = results.stream().filter(r -> r.source() == src && r.paths() == n).toList();
                out.printf(Locale.ROOT, "   N = %d%n", n);
                for (int c = 0; c < configs.size(); c++) {
                    final int ci = c;
                    double[] loss = rs.stream().mapToDouble(r -> r.loss()[ci]).toArray();
                    double[] diff = rs.stream().mapToDouble(r -> r.loss()[ci] - r.loss()[baseline]).toArray();
                    double better = rs.stream().filter(r -> r.loss()[ci] < r.loss()[baseline]).count() / (double) rs.size();
                    double[] bias = rs.stream().mapToDouble(r -> r.inSampleBias()[ci]).toArray();
                    double[] rms = rs.stream().mapToDouble(r -> r.continuationRms()[ci]).toArray();
                    double absB = 0, sumB = 0, cErr2 = 0;
                    int nb = 0, nonThreshold = 0, nc = 0;
                    for (Replication r : rs) {
                        for (int j = 0; j < dates; j++) {
                            double e = r.boundaryError()[ci][j];
                            if (!Double.isNaN(e)) {
                                absB += Math.abs(e);
                                sumB += e;
                                nb++;
                            }
                            if (r.switchCount()[ci][j] != 1) nonThreshold++;
                            double ce = r.continuationErrorAtBoundary()[ci][j];
                            if (!Double.isNaN(ce)) {
                                cErr2 += ce * ce;
                                nc++;
                            }
                        }
                    }
                    double nonThresholdShare = nonThreshold / (double) (rs.size() * dates);
                    out.printf(Locale.ROOT, "     [%d] %-10s %.3e (%.1e) | %.2e | %.2e   diff %s   |r^-r*| %.2e%s   C-err %.2e   in-sample %+.2e (%.1e)%n",
                            c, configs.get(c).name(), mean(loss), se(loss), quantile(loss, 0.5), quantile(loss, 0.9),
                            c == baseline ? "        (baseline)        "
                                    : String.format(Locale.ROOT, "%+.2e (%.1e, %3.0f%%)", mean(diff), se(diff), 100 * better),
                            nb == 0 ? Double.NaN : absB / nb,
                            nonThresholdShare > 0 ? String.format(Locale.ROOT, " (%.0f%% non-threshold)", 100 * nonThresholdShare) : "",
                            nc == 0 ? Double.NaN : Math.sqrt(cErr2 / nc), mean(bias), se(bias));
                    summary.row(id, src, n, configs.get(c).name(), mean(loss), se(loss), quantile(loss, 0.5),
                            quantile(loss, 0.9), c == baseline ? "" : mean(diff), c == baseline ? "" : se(diff),
                            c == baseline ? "" : better, nb == 0 ? "" : absB / nb, nb == 0 ? "" : sumB / nb,
                            nonThresholdShare, nc == 0 ? "" : Math.sqrt(cErr2 / nc), mean(rms), mean(bias), se(bias));
                }
            }
        }

        if (cfg.writeCsv) {
            java.nio.file.Path dir = cfg.outputDir.resolve(id);
            raw.write(dir);
            summary.write(dir);
            Map<String, String> extra = new LinkedHashMap<>();
            extra.put("kind", "method study (cluster vs LSM vs exact solution)");
            extra.put("v_star", String.valueOf(vStar));
            extra.put("r_star", Arrays.toString(rStar));
            for (int c = 0; c < configs.size(); c++) extra.put("config." + configs.get(c).name(), configs.get(c).settings().label());
            RunMetadata.write(dir, id, cfg.entries(), extra);
            out.println("CSV written to " + dir.toAbsolutePath());
        }
    }

    private static Replication replicate(PathSource src, int n, int rep, long seed, ExactPathGenerator exact, Setup euler,
                                         List<Config> configs, CirQuadratureBermudan quad, ExperimentConfig cfg,
                                         CirQuadratureBermudan.Settings qs, double vStar,
                                         CirQuadratureBermudan.Result optimal, double[] rStar) {
        PathDataDouble train = src == PathSource.EXACT ? exact.simulate(seed, n, 1)
                : euler.simulateDouble(euler.source(seed), n, 1);
        int k = configs.size();
        int dates = rStar.length;
        double[] loss = new double[k];
        double[] inSample = new double[k];
        double[] rms = new double[k];
        double[][] bErr = new double[k][dates];
        int[][] switches = new int[k][dates];
        double[][] cErr = new double[k][dates];
        for (int c = 0; c < k; c++) {
            InductionSettings s = configs.get(c).settings();
            InductionResultDouble ind = BackwardInductionDouble.run(train, TrainingStructure.fromDouble(train, s), s);
            ExercisePolicy policy = policyOf(ind);
            CirQuadratureBermudan.Result v = quad.valuePolicy(cfg.option, policy, qs);
            loss[c] = vStar - v.price();
            inSample[c] = ind.price.mean() - vStar;
            double sq = 0.0;
            int cnt = 0;
            for (int j = 0; j < dates; j++) {
                double[] sw = v.exerciseBoundaries()[j];
                switches[c][j] = sw.length;
                bErr[c][j] = Double.NaN;
                double best = Double.POSITIVE_INFINITY;
                for (double x : sw) {
                    if (Math.abs(x - rStar[j]) < Math.abs(best)) best = x - rStar[j];
                }
                if (sw.length > 0) bErr[c][j] = best;
                cErr[c][j] = estimate(ind.rules[j], rStar[j], train, j) - optimal.continuationValue(j, rStar[j]);
                for (int p = 0; p < train.paths(); p++) {
                    double r = train.rate[j][p];
                    if (r > qs.rateMax()) continue;
                    double e = estimate(ind.rules[j], r, train, j) - optimal.continuationValue(j, r);
                    if (!Double.isNaN(e)) {
                        sq += e * e;
                        cnt++;
                    }
                }
            }
            rms[c] = cnt == 0 ? Double.NaN : Math.sqrt(sq / cnt);
        }
        return new Replication(src, n, rep, loss, bErr, switches, cErr, rms, inSample);
    }

    /** Exercise rule of a trained induction: exercise iff {@code E > 0} and {@code E > C^(state)}; last date {@code C^ = 0}. */
    static ExercisePolicy policyOf(InductionResultDouble ind) {
        return (j, rate, bondValue, e) -> {
            if (!(e > 0.0)) return false;
            DateRuleDouble rule = ind.rules[j];
            if (rule == null) return true;
            double c = continuation(rule, rate, bondValue);
            return e > c;
        };
    }

    static double continuation(DateRuleDouble rule, double rate, double bondValue) {
        if (rule instanceof ClusterRuleDouble cr) {
            return cr.clusterCount() == 0 ? Double.NaN : cr.mean()[cr.assign(rate)];
        }
        return ((RegressionRuleDouble) rule).predict(rate, bondValue);
    }

    private static double estimate(DateRuleDouble rule, double rate, PathDataDouble d, int j) {
        // bond value at an arbitrary rate is only needed for the bond-value regression variable; recompute exactly
        double bond = Double.NaN;
        if (rule instanceof RegressionRuleDouble rr && rr.variable() != lab.cir.induction.RegressionVariable.SHORT_RATE) {
            throw new UnsupportedOperationException("method study supports the short-rate regression variable only");
        }
        return continuation(rule, rate, bond);
    }

    static List<Config> configs(ExperimentConfig cfg) {
        List<Config> list = new ArrayList<>();
        InductionSettings base = cfg.inductionSettings(ContinuationMethod.QUANTILE_CLUSTER)
                .withPathSelection(PathSelection.ALL);
        for (String k : EuropeanValidation.getOr(cfg, "methodstudy.clusters", "5,10,20,40").split(",")) {
            int kk = Integer.parseInt(k.trim());
            list.add(new Config("CL" + kk, base.withClusterCount(kk)));
        }
        for (String spec : EuropeanValidation.getOr(cfg, "methodstudy.lsm", "2:ALL,2:IN_THE_MONEY,3:IN_THE_MONEY,4:IN_THE_MONEY")
                .split(",")) {
            String[] parts = spec.trim().split(":");
            int degree = Integer.parseInt(parts[0]);
            PathSelection sel = PathSelection.valueOf(parts[1]);
            list.add(new Config("LSM" + degree + "-" + (sel == PathSelection.ALL ? "ALL" : "ITM"),
                    base.withMethod(ContinuationMethod.LEAST_SQUARES).withPathSelection(sel).withBasisDegree(degree)));
        }
        return list;
    }

    private static int baselineIndex(List<Config> configs, String name) {
        for (int i = 0; i < configs.size(); i++) if (configs.get(i).name().equals(name)) return i;
        return 0;
    }

    private static int[] ints(String csv) {
        return Arrays.stream(csv.split(",")).mapToInt(s -> Integer.parseInt(s.trim())).toArray();
    }

    static double mean(double[] x) {
        double s = 0.0;
        for (double v : x) s += v;
        return s / x.length;
    }

    static double se(double[] x) {
        if (x.length < 2) return Double.NaN;
        double m = mean(x);
        double ss = 0.0;
        for (double v : x) ss += (v - m) * (v - m);
        return Math.sqrt(ss / (x.length - 1) / x.length);
    }

    static double quantile(double[] x, double q) {
        double[] s = x.clone();
        Arrays.sort(s);
        double pos = q * (s.length - 1);
        int lo = (int) Math.floor(pos);
        int hi = Math.min(s.length - 1, lo + 1);
        return s[lo] + (pos - lo) * (s[hi] - s[lo]);
    }
}

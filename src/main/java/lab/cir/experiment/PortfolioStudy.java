package lab.cir.experiment;

import lab.cir.analytic.CirZeroBondDecimal;
import lab.cir.analytic.CirZeroBondDouble;
import lab.cir.analytic.RemainingBondDecimal;
import lab.cir.analytic.RemainingBondDouble;
import lab.cir.config.ExperimentConfig;
import lab.cir.induction.BackwardInductionDecimal;
import lab.cir.induction.BackwardInductionDouble;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionResultDecimal;
import lab.cir.induction.InductionResultDouble;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.TrainingStructure;
import lab.cir.model.BermudanCallOption;
import lab.cir.numeric.DecimalMath;
import lab.cir.random.NormalSource;
import lab.cir.report.CsvTable;
import lab.cir.report.Fmt;
import lab.cir.report.Measured;
import lab.cir.report.RunMetadata;
import lab.cir.simulation.PathDataDecimal;
import lab.cir.simulation.PathDataDouble;
import lab.cir.simulation.PathGeneratorDecimal;
import lab.cir.simulation.PathGeneratorDouble;
import lab.cir.simulation.TimeGrid;
import lab.cir.stats.SampleStatsDecimal;
import lab.cir.stats.SampleStatsDouble;

import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.LinkedHashMap;

/**
 * Portfolio of Bermudan (and European) calls on the configured bond, heterogeneous in strike, notional scaling and
 * exercise dates. Two strategies, both in-sample (training paths = pricing paths) with quantile clusters:
 * <ol>
 *   <li>independent valuation: every contract simulates its own paths (seed {@code seed + 7919 c}); the portfolio
 *       standard error is {@code sqrt(sum SE_c^2)} (independent estimates);</li>
 *   <li>shared market simulation: one set of CIR paths on the union of exercise dates; exercise values, clusters and
 *       policies remain per contract; the portfolio standard error comes from the path-wise aggregated cashflows.</li>
 * </ol>
 * No result caching between contracts: every contract runs its own exercise values and induction.
 */
public final class PortfolioStudy {

    private PortfolioStudy() {
    }

    public static List<BermudanCallOption> portfolio(BermudanCallOption template, int contracts) {
        String[] strikes = {"92", "96", "100", "104", "108"};
        String[] scales = {"1", "2.5", "0.5", "4", "1.5"};
        int m = template.exerciseDateCount();
        List<BermudanCallOption> out = new ArrayList<>();
        for (int c = 0; c < contracts; c++) {
            BigDecimal scale = new BigDecimal(scales[(c / 5) % scales.length]);
            BigDecimal strike = new BigDecimal(strikes[c % strikes.length]);
            // exercise date subsets: all, drop first, drop last, single (European) at first, alternating
            List<BigDecimal> dates = new ArrayList<>();
            int pattern = (c / 25) % 5;
            for (int j = 0; j < m; j++) {
                boolean keep = switch (pattern) {
                    case 0 -> true;
                    case 1 -> j > 0;
                    case 2 -> j < m - 1;
                    case 3 -> j == 0;
                    default -> j % 2 == 0;
                };
                if (keep) dates.add(template.exerciseDate(j));
            }
            out.add(new BermudanCallOption(template.underlying().scaled(scale), strike.multiply(scale), dates));
        }
        return out;
    }

    public static void run(ExperimentConfig cfg, int contracts, PrintStream out) throws IOException {
        String id = RunMetadata.newExperimentId("portfolio-" + contracts, cfg.label);
        List<BermudanCallOption> book = portfolio(cfg.option, contracts);
        InductionSettings is = cfg.inductionSettings(ContinuationMethod.QUANTILE_CLUSTER);
        DecimalMath math = new DecimalMath(cfg.precision);
        MathContext mc = math.mc();
        boolean withDecimal = !"false".equalsIgnoreCase(EuropeanValidation.getOr(cfg, "portfolio.decimal", "true"));
        int paths = cfg.trainingPaths;
        out.println("== Portfolio study: " + contracts + " contracts, " + paths + " paths, " + cfg.stepsPerYear
                + " steps/year, clusters K=" + cfg.clusterCount + ", " + cfg.precision.precision() + " digits  id=" + id);

        CsvTable t = new CsvTable("portfolio", "experiment_id", "contracts", "paths", "strategy", "variant", "portfolio_value",
                "portfolio_se", "total_ms", "ms_per_contract", "contracts_per_second", "allocated_bytes", "peak_heap_used_bytes");
        Map<String, Object[]> results = new LinkedHashMap<>();

        // ---------------- independent valuation
        Measured<double[]> indD = Measured.of(() -> {
            double sum = 0, var = 0;
            for (int c = 0; c < book.size(); c++) {
                Setup s = new Setup(cfg.model, book.get(c), cfg.stepsPerYear, cfg.normalsPerStep, cfg.precision, cfg.zeroBondForm);
                PathDataDouble d = s.simulateDouble(s.source(cfg.trainingSeed + 7919L * c), paths, cfg.threads);
                InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, is), is);
                sum += r.price.mean();
                var += r.price.standardError() * r.price.standardError();
            }
            return new double[]{sum, Math.sqrt(var)};
        });
        record(out, t, id, contracts, paths, "independent", "double", indD.value()[0], indD.value()[1], indD, results);
        if (withDecimal) {
            resetPeak();
            Measured<BigDecimal[]> indB = Measured.of(() -> {
                BigDecimal sum = BigDecimal.ZERO, var = BigDecimal.ZERO;
                for (int c = 0; c < book.size(); c++) {
                    Setup s = new Setup(cfg.model, book.get(c), cfg.stepsPerYear, cfg.normalsPerStep, cfg.precision, cfg.zeroBondForm);
                    PathDataDecimal d = s.simulateDecimal(s.source(cfg.trainingSeed + 7919L * c), paths, cfg.threads);
                    InductionResultDecimal r = BackwardInductionDecimal.run(d, TrainingStructure.fromDecimal(d, is), is, s.math);
                    sum = sum.add(r.price.mean(), mc);
                    var = var.add(r.price.standardError().multiply(r.price.standardError(), mc), mc);
                }
                return new BigDecimal[]{sum, math.sqrt(var)};
            });
            record(out, t, id, contracts, paths, "independent", "decimal", indB.value()[0].doubleValue(),
                    indB.value()[1].doubleValue(), indB, results);
            results.put("independent-decimal-exact", indB.value());
        }

        // ---------------- shared market simulation
        TreeSet<BigDecimal> union = new TreeSet<>();
        for (BermudanCallOption o : book) for (int j = 0; j < o.exerciseDateCount(); j++) union.add(o.exerciseDate(j));
        List<BigDecimal> unionDates = new ArrayList<>(union);
        TimeGrid grid = new TimeGrid(unionDates, cfg.stepsPerYear, cfg.normalsPerStep);
        NormalSource src = new NormalSource(cfg.trainingSeed, grid.normalsPerPath());
        resetPeak();
        Measured<double[]> shD = Measured.of(() -> {
            PathDataDouble market = new PathGeneratorDouble(cfg.model, grid).simulate(src, paths, cfg.threads);
            CirZeroBondDouble zb = new CirZeroBondDouble(cfg.model, cfg.zeroBondForm);
            double[] pathSum = new double[paths];
            double sum = 0;
            for (BermudanCallOption o : book) {
                PathDataDouble d = new PathDataDouble(market, indices(unionDates, o));
                RemainingBondDouble bond = new RemainingBondDouble(o, zb);
                fillExerciseValues(d, bond);
                InductionResultDouble r = BackwardInductionDouble.run(d, TrainingStructure.fromDouble(d, is), is);
                sum += r.price.mean();
                for (int p = 0; p < paths; p++) pathSum[p] += r.pathValue[p];
            }
            return new double[]{sum, SampleStatsDouble.of(pathSum).standardError()};
        });
        record(out, t, id, contracts, paths, "shared", "double", shD.value()[0], shD.value()[1], shD, results);
        if (withDecimal) {
            resetPeak();
            Measured<BigDecimal[]> shB = Measured.of(() -> {
                PathDataDecimal market = new PathGeneratorDecimal(cfg.model, grid, math).simulate(src, paths, cfg.threads);
                CirZeroBondDecimal zb = new CirZeroBondDecimal(cfg.model, math, cfg.zeroBondForm);
                BigDecimal[] pathSum = new BigDecimal[paths];
                java.util.Arrays.fill(pathSum, BigDecimal.ZERO);
                BigDecimal sum = BigDecimal.ZERO;
                for (BermudanCallOption o : book) {
                    PathDataDecimal d = new PathDataDecimal(market, indices(unionDates, o));
                    RemainingBondDecimal bond = new RemainingBondDecimal(o, zb);
                    fillExerciseValues(d, bond);
                    InductionResultDecimal r = BackwardInductionDecimal.run(d, TrainingStructure.fromDecimal(d, is), is, math);
                    sum = sum.add(r.price.mean(), mc);
                    for (int p = 0; p < paths; p++) pathSum[p] = pathSum[p].add(r.pathValue[p], mc);
                }
                return new BigDecimal[]{sum, SampleStatsDecimal.of(pathSum, math).standardError()};
            });
            record(out, t, id, contracts, paths, "shared", "decimal", shB.value()[0].doubleValue(), shB.value()[1].doubleValue(),
                    shB, results);
            results.put("shared-decimal-exact", shB.value());
            BigDecimal dIndep = ((BigDecimal[]) results.get("independent-decimal-exact"))[0].subtract(new BigDecimal(indD.value()[0]), mc);
            BigDecimal dShared = shB.value()[0].subtract(new BigDecimal(shD.value()[0]), mc);
            out.printf("   portfolio Delta (dec - dbl): independent %s, shared %s;  shared - independent (double) %s%n",
                    Fmt.de(dIndep), Fmt.de(dShared), Fmt.e(shD.value()[0] - indD.value()[0]));
        }
        if (cfg.writeCsv) {
            java.nio.file.Path dir = cfg.outputDir.resolve(id);
            t.write(dir);
            Map<String, String> extra = new LinkedHashMap<>();
            extra.put("kind", "portfolio study (in-sample, quantile clusters)");
            extra.put("contracts", String.valueOf(contracts));
            RunMetadata.write(dir, id, cfg.entries(), extra);
            out.println("CSV written to " + dir.toAbsolutePath());
        }
    }

    private static int[] indices(List<BigDecimal> union, BermudanCallOption o) {
        int[] idx = new int[o.exerciseDateCount()];
        for (int j = 0; j < idx.length; j++) idx[j] = union.indexOf(o.exerciseDate(j));
        return idx;
    }

    private static void fillExerciseValues(PathDataDouble d, RemainingBondDouble bond) {
        for (int j = 0; j < d.dates(); j++) {
            for (int p = 0; p < d.paths(); p++) {
                double b = bond.bondValue(j, d.rate[j][p]);
                d.bondValue[j][p] = b;
                d.exerciseValue[j][p] = bond.exerciseValue(b);
            }
        }
    }

    private static void fillExerciseValues(PathDataDecimal d, RemainingBondDecimal bond) {
        for (int j = 0; j < d.dates(); j++) {
            for (int p = 0; p < d.paths(); p++) {
                BigDecimal b = bond.bondValue(j, d.rate[j][p]);
                d.bondValue[j][p] = b;
                d.exerciseValue[j][p] = bond.exerciseValue(b);
            }
        }
    }

    private static void resetPeak() {
        System.gc();
        for (java.lang.management.MemoryPoolMXBean b : java.lang.management.ManagementFactory.getMemoryPoolMXBeans()) {
            b.resetPeakUsage();
        }
    }

    private static long peakHeap() {
        long s = 0;
        for (java.lang.management.MemoryPoolMXBean b : java.lang.management.ManagementFactory.getMemoryPoolMXBeans()) {
            if (b.getType() == java.lang.management.MemoryType.HEAP && b.getPeakUsage() != null) s += b.getPeakUsage().getUsed();
        }
        return s;
    }

    private static void record(PrintStream out, CsvTable t, String id, int contracts, int paths, String strategy, String variant,
                               double value, double se, Measured<?> m, Map<String, Object[]> results) {
        double ms = m.millis();
        long peak = peakHeap();
        out.printf("   %-12s %-8s value %16s  SE %10s  total %10.1f ms  %8.2f ms/contract  %9.1f contracts/s  alloc %9s MB  peak heap %8s MB%n",
                strategy, variant, Fmt.g(value), Fmt.e(se), ms, ms / contracts, contracts / (ms / 1000.0),
                Fmt.mb(m.allocatedBytes()), Fmt.mb(peak));
        t.row(id, contracts, paths, strategy, variant, value, se, ms, ms / contracts, contracts / (ms / 1000.0),
                m.allocatedBytes(), peak);
    }
}

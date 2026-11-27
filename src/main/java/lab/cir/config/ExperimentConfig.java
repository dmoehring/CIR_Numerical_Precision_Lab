package lab.cir.config;

import lab.cir.analytic.ZeroBondForm;
import lab.cir.induction.ContinuationMethod;
import lab.cir.induction.InductionSettings;
import lab.cir.induction.PathSelection;
import lab.cir.induction.RegressionVariable;
import lab.cir.induction.SummationMode;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;
import lab.cir.model.CirParameters;
import lab.cir.model.FellerMode;
import lab.cir.numeric.PrecisionSettings;

import java.io.IOException;
import java.io.Reader;
import java.math.RoundingMode;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.TreeMap;

/**
 * Experiment configuration read from a properties file (plus {@code key=value} overrides). All business inputs
 * are decimal strings; both numerical variants are derived from the same strings.
 */
public final class ExperimentConfig {

    private final Properties raw;

    public final CirParameters model;
    public final BermudanCallOption option;
    public final int stepsPerYear;
    public final int normalsPerStep;
    public final int trainingPaths;
    public final long trainingSeed;
    public final int valuationPaths;
    public final long valuationSeed;
    public final int threads;
    public final int clusterCount;
    public final int minPathsPerCluster;
    public final PathSelection pathSelection;
    public final RegressionVariable regressionVariable;
    public final int basisDegree;
    public final List<ContinuationMethod> methods;
    public final PrecisionSettings precision;
    public final ZeroBondForm zeroBondForm;
    public final int tracePaths;
    public final boolean kahanControl;
    public final Path outputDir;
    public final boolean writeCsv;
    public final String label;

    private ExperimentConfig(Properties p) {
        this.raw = p;
        FellerMode mode = FellerMode.valueOf(get(p, "cir.fellerMode", "STRICT"));
        this.model = CirParameters.of(req(p, "cir.meanReversionSpeed"), req(p, "cir.longTermMean"),
                req(p, "cir.volatility"), req(p, "cir.initialRate"), mode);
        CashflowSchedule bond = p.containsKey("bond.cashflowTimes")
                ? new CashflowSchedule(decimals(req(p, "bond.cashflowTimes")), decimals(req(p, "bond.cashflowAmounts")))
                : CashflowSchedule.couponBond(req(p, "bond.notional"), req(p, "bond.couponRate"),
                req(p, "bond.maturity"), Integer.parseInt(get(p, "bond.couponFrequency", "1")));
        this.option = new BermudanCallOption(bond, new java.math.BigDecimal(req(p, "option.strike").trim()),
                decimals(req(p, "option.exerciseDates")));
        this.stepsPerYear = Integer.parseInt(req(p, "simulation.stepsPerYear"));
        this.normalsPerStep = Integer.parseInt(get(p, "simulation.normalsPerStep", "1"));
        this.trainingPaths = Integer.parseInt(req(p, "simulation.paths"));
        this.trainingSeed = Long.parseLong(req(p, "simulation.seed"));
        this.valuationPaths = Integer.parseInt(get(p, "valuation.paths", String.valueOf(trainingPaths)));
        this.valuationSeed = Long.parseLong(get(p, "valuation.seed", String.valueOf(trainingSeed + 1_000_003L)));
        this.threads = Integer.parseInt(get(p, "simulation.threads", "1"));
        this.clusterCount = Integer.parseInt(get(p, "cluster.count", "10"));
        this.minPathsPerCluster = Integer.parseInt(get(p, "cluster.minPathsPerCluster", "20"));
        this.pathSelection = PathSelection.valueOf(get(p, "continuation.pathSelection", "ALL"));
        this.regressionVariable = RegressionVariable.valueOf(get(p, "lsm.basisVariable", "SHORT_RATE"));
        this.basisDegree = Integer.parseInt(get(p, "lsm.basisDegree", "2"));
        this.methods = Arrays.stream(get(p, "continuation.methods", "QUANTILE_CLUSTER,LEAST_SQUARES").split(","))
                .map(String::trim).filter(s -> !s.isEmpty()).map(ContinuationMethod::valueOf).toList();
        this.precision = new PrecisionSettings(Integer.parseInt(get(p, "decimal.precision", "34")),
                RoundingMode.valueOf(get(p, "decimal.roundingMode", "HALF_EVEN")),
                Integer.parseInt(get(p, "decimal.guardDigits", "10")));
        this.zeroBondForm = ZeroBondForm.valueOf(get(p, "zerobond.form", "LOG_DIFFERENCE"));
        this.tracePaths = Integer.parseInt(get(p, "diagnostics.tracePaths", "200"));
        this.kahanControl = Boolean.parseBoolean(get(p, "control.kahan", "true"));
        this.outputDir = Path.of(get(p, "output.dir", "results"));
        this.writeCsv = Boolean.parseBoolean(get(p, "output.csv", "true"));
        this.label = get(p, "experiment.label", "experiment");
        if (trainingPaths < 2 || valuationPaths < 2) {
            throw new IllegalArgumentException("need at least 2 training and valuation paths");
        }
    }

    public static ExperimentConfig load(Path file, Map<String, String> overrides) throws IOException {
        Properties p = new Properties();
        try (Reader r = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            p.load(r);
        }
        overrides.forEach(p::setProperty);
        return new ExperimentConfig(p);
    }

    public static ExperimentConfig of(Properties p) {
        Properties copy = new Properties();
        copy.putAll(p);
        return new ExperimentConfig(copy);
    }

    /** Copy with one key replaced (single-factor sensitivity series). */
    public ExperimentConfig with(String key, String value) {
        Properties copy = new Properties();
        copy.putAll(raw);
        copy.setProperty(key, value);
        return new ExperimentConfig(copy);
    }

    public InductionSettings inductionSettings(ContinuationMethod method) {
        return new InductionSettings(method, pathSelection, clusterCount, minPathsPerCluster, regressionVariable,
                SummationMode.NAIVE, basisDegree);
    }

    /** Sorted view of all configuration entries (for metadata files). */
    public Map<String, String> entries() {
        Map<String, String> m = new TreeMap<>();
        for (String k : raw.stringPropertyNames()) m.put(k, raw.getProperty(k));
        return m;
    }

    public String get(String key) {
        return raw.getProperty(key);
    }

    private static String req(Properties p, String key) {
        String v = p.getProperty(key);
        if (v == null || v.isBlank()) {
            throw new IllegalArgumentException("missing configuration key " + key);
        }
        return v.trim();
    }

    private static String get(Properties p, String key, String def) {
        String v = p.getProperty(key);
        return v == null || v.isBlank() ? def : v.trim();
    }

    private static List<java.math.BigDecimal> decimals(String csv) {
        List<java.math.BigDecimal> out = new ArrayList<>();
        for (String s : csv.split(",")) {
            if (!s.isBlank()) out.add(new java.math.BigDecimal(s.trim()));
        }
        return out;
    }
}

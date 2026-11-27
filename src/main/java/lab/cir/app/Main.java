package lab.cir.app;

import lab.cir.config.ExperimentConfig;
import lab.cir.experiment.ComparisonExperiment;
import lab.cir.experiment.EuropeanValidation;
import lab.cir.experiment.FdConvergenceStudy;
import lab.cir.experiment.MethodStudy;
import lab.cir.experiment.PortfolioStudy;
import lab.cir.experiment.ReferenceStudy;
import lab.cir.experiment.SensitivityStudy;
import lab.cir.report.ComparisonReport;
import lab.cir.report.RunMetadata;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Command line entry point.
 *
 * <pre>
 *   compare     &lt;config&gt; [key=value ...]                       double vs BigDecimal, Mode A/B, cluster + LSM
 *   european    &lt;config&gt; [key=value ...]                       European control cases vs analytic references
 *   reference   &lt;config&gt; [key=value ...]                       quadrature/FD Bermudan reference + Euler/time-step series
 *   fdstudy     &lt;config&gt; [key=value ...]                       accuracy of the quadrature and FD references
 *   methodstudy &lt;config&gt; [key=value ...]                       cluster vs LSM policies valued exactly (seeds x paths)
 *   sensitivity &lt;config&gt; &lt;key&gt; &lt;v1,v2,...&gt; [key=value ...]   one-factor series (paths, steps, clusters, precision ...)
 *   portfolio   &lt;config&gt; &lt;contracts&gt; [key=value ...]          independent vs shared simulation
 * </pre>
 */
public final class Main {

    private Main() {
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            usage();
            return;
        }
        String command = args[0];
        Path configFile = Path.of(args[1]);
        int fixedArgs = switch (command) {
            case "sensitivity" -> 4;
            case "portfolio" -> 3;
            default -> 2;
        };
        if (args.length < fixedArgs) {
            usage();
            return;
        }
        Map<String, String> overrides = overrides(Arrays.asList(args).subList(fixedArgs, args.length));
        ExperimentConfig cfg = ExperimentConfig.load(configFile, overrides);
        switch (command) {
            case "compare" -> {
                String id = RunMetadata.newExperimentId("compare", cfg.label);
                ComparisonExperiment.Run run = ComparisonExperiment.run(cfg, true);
                ComparisonReport.print(run, id, System.out);
                if (cfg.writeCsv) {
                    Path dir = cfg.outputDir.resolve(id);
                    ComparisonReport.writeCsv(run, id, dir);
                    System.out.println("CSV written to " + dir.toAbsolutePath());
                }
            }
            case "european" -> EuropeanValidation.run(cfg, System.out);
            case "reference" -> ReferenceStudy.run(cfg, System.out);
            case "fdstudy" -> FdConvergenceStudy.run(cfg, System.out);
            case "methodstudy" -> MethodStudy.run(cfg, System.out);
            case "sensitivity" -> SensitivityStudy.run(cfg, args[2], List.of(args[3].split(",")), System.out);
            case "portfolio" -> PortfolioStudy.run(cfg, Integer.parseInt(args[2]), System.out);
            default -> usage();
        }
    }

    private static Map<String, String> overrides(List<String> args) {
        Map<String, String> m = new LinkedHashMap<>();
        for (String a : args) {
            int i = a.indexOf('=');
            if (i <= 0) throw new IllegalArgumentException("override must be key=value: " + a);
            m.put(a.substring(0, i).trim(), a.substring(i + 1).trim());
        }
        return m;
    }

    private static void usage() {
        System.out.println("""
                usage: java -jar target/cir-numerical-precision-lab.jar <command> <config.properties> [args] [key=value ...]
                  compare     <config>                          full double vs BigDecimal comparison (console + CSV)
                  european    <config>                          European control cases vs analytic references
                  reference   <config>                          quadrature/FD Bermudan reference, time-step series
                  fdstudy     <config>                          accuracy of quadrature and finite-difference references
                  methodstudy <config>                          cluster vs LSM: trained policies valued exactly, many seeds
                  sensitivity <config> <key> <v1,v2,...>        one-factor series, e.g. simulation.paths 1000,4000,16000
                  portfolio   <config> <contracts>              independent valuation vs shared market simulation
                JMH: java -cp target/cir-numerical-precision-lab.jar org.openjdk.jmh.Main lab.cir.bench
                """);
    }
}

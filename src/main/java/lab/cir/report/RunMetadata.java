package lab.cir.report;

import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

/** Experiment id and runtime environment (Java, JVM args, CPU, OS, heap, GC, library versions, git commit). */
public final class RunMetadata {

    private RunMetadata() {
    }

    public static String newExperimentId(String kind, String label) {
        String ts = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd-HHmmss-SSS"));
        return ts + "-" + kind + (label == null || label.isBlank() ? "" : "-" + label.replaceAll("[^A-Za-z0-9_.-]", "_"));
    }

    public static Map<String, String> environment() {
        Map<String, String> m = new LinkedHashMap<>();
        m.put("java.version", System.getProperty("java.version"));
        m.put("java.vm.name", System.getProperty("java.vm.name") + " " + System.getProperty("java.vm.version"));
        m.put("jvm.args", String.join(" ", ManagementFactory.getRuntimeMXBean().getInputArguments()));
        m.put("os", System.getProperty("os.name") + " " + System.getProperty("os.version") + " " + System.getProperty("os.arch"));
        m.put("cpu", cpuModel());
        m.put("cpu.available", String.valueOf(Runtime.getRuntime().availableProcessors()));
        m.put("heap.max.mb", String.valueOf(Runtime.getRuntime().maxMemory() / (1024 * 1024)));
        StringBuilder gc = new StringBuilder();
        for (GarbageCollectorMXBean b : ManagementFactory.getGarbageCollectorMXBeans()) {
            if (!gc.isEmpty()) gc.append("; ");
            gc.append(b.getName());
        }
        m.put("gc", gc.toString());
        Properties versions = new Properties();
        try (InputStream in = RunMetadata.class.getResourceAsStream("/lab-versions.properties")) {
            if (in != null) versions.load(in);
        } catch (IOException ignored) {
            // versions stay empty
        }
        versions.stringPropertyNames().stream().sorted().forEach(k -> m.put("lib." + k, versions.getProperty(k)));
        m.put("git.commit", gitCommit());
        return m;
    }

    public static void write(Path dir, String experimentId, Map<String, String> config, Map<String, String> extra) throws IOException {
        Files.createDirectories(dir);
        try (Writer w = Files.newBufferedWriter(dir.resolve("run_metadata.properties"), StandardCharsets.UTF_8)) {
            w.write("experiment.id=" + experimentId + "\n");
            for (Map.Entry<String, String> e : environment().entrySet()) {
                w.write("env." + e.getKey() + "=" + e.getValue() + "\n");
            }
            for (Map.Entry<String, String> e : extra.entrySet()) {
                w.write("run." + e.getKey() + "=" + e.getValue() + "\n");
            }
            for (Map.Entry<String, String> e : config.entrySet()) {
                w.write("config." + e.getKey() + "=" + e.getValue() + "\n");
            }
        }
    }

    private static String cpuModel() {
        try {
            Path p = Path.of("/proc/cpuinfo");
            if (Files.exists(p)) {
                for (String line : Files.readAllLines(p)) {
                    if (line.startsWith("model name")) {
                        return line.substring(line.indexOf(':') + 1).trim();
                    }
                }
            }
        } catch (IOException ignored) {
            // fall through
        }
        return System.getProperty("os.arch");
    }

    private static String gitCommit() {
        try {
            Process p = new ProcessBuilder("git", "rev-parse", "--short", "HEAD").redirectErrorStream(true).start();
            String head = new String(p.getInputStream().readAllBytes(), StandardCharsets.UTF_8).trim();
            if (!p.waitFor(5, TimeUnit.SECONDS) || p.exitValue() != 0) return "unavailable";
            Process s = new ProcessBuilder("git", "status", "--porcelain", "--untracked-files=no").redirectErrorStream(true).start();
            String dirty = new String(s.getInputStream().readAllBytes(), StandardCharsets.UTF_8).trim();
            s.waitFor(5, TimeUnit.SECONDS);
            return head + (dirty.isEmpty() ? "" : "-dirty");
        } catch (IOException e) {
            return "unavailable";
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return "unavailable";
        }
    }
}

package lab.cir.report;

import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Minimal CSV table. BigDecimal values are written with {@link BigDecimal#toString()} (full precision, no prior
 * double conversion); doubles with {@link Double#toString(double)} (shortest round-trip representation).
 */
public final class CsvTable {

    private final String name;
    private final List<String> header;
    private final List<List<String>> rows = new ArrayList<>();

    public CsvTable(String name, String... header) {
        this.name = name;
        this.header = List.of(header);
    }

    public CsvTable row(Object... values) {
        if (values.length != header.size()) {
            throw new IllegalArgumentException(name + ": expected " + header.size() + " values, got " + values.length);
        }
        List<String> r = new ArrayList<>(values.length);
        for (Object v : values) r.add(format(v));
        rows.add(r);
        return this;
    }

    public String name() { return name; }
    public int size() { return rows.size(); }

    public void write(Path dir) throws IOException {
        Files.createDirectories(dir);
        try (Writer w = Files.newBufferedWriter(dir.resolve(name + ".csv"), StandardCharsets.UTF_8)) {
            w.write(String.join(",", header));
            w.write('\n');
            for (List<String> r : rows) {
                w.write(String.join(",", r));
                w.write('\n');
            }
        }
    }

    static String format(Object v) {
        if (v == null) return "";
        if (v instanceof BigDecimal b) return b.toString();
        if (v instanceof Double d) return d.isNaN() ? "NaN" : Double.toString(d);
        String s = v.toString();
        if (s.contains(",") || s.contains("\"") || s.contains("\n")) {
            return '"' + s.replace("\"", "\"\"") + '"';
        }
        return s;
    }
}

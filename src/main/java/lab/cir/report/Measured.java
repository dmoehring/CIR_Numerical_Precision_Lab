package lab.cir.report;

import java.lang.management.ManagementFactory;
import java.util.function.Supplier;

/**
 * Orientation timing with {@link System#nanoTime()} and bytes allocated by the calling thread
 * ({@code com.sun.management.ThreadMXBean}). Not a substitute for JMH: no warm-up, single run.
 * Allocation counts only the calling thread, i.e. they are complete only for {@code simulation.threads=1}.
 */
public record Measured<T>(T value, long nanos, long allocatedBytes) {

    private static final com.sun.management.ThreadMXBean THREADS =
            ManagementFactory.getThreadMXBean() instanceof com.sun.management.ThreadMXBean t ? t : null;

    static {
        if (THREADS != null && THREADS.isThreadAllocatedMemorySupported()) {
            THREADS.setThreadAllocatedMemoryEnabled(true);
        }
    }

    public static <T> Measured<T> of(Supplier<T> task) {
        long a0 = allocated();
        long t0 = System.nanoTime();
        T v = task.get();
        long t1 = System.nanoTime();
        long a1 = allocated();
        return new Measured<>(v, t1 - t0, a1 < 0 || a0 < 0 ? -1 : a1 - a0);
    }

    public static Measured<Void> run(Runnable task) {
        return of(() -> {
            task.run();
            return null;
        });
    }

    public double millis() {
        return nanos / 1e6;
    }

    private static long allocated() {
        return THREADS == null ? -1 : THREADS.getCurrentThreadAllocatedBytes();
    }
}

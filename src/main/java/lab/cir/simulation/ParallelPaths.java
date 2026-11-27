package lab.cir.simulation;

import java.util.concurrent.ForkJoinPool;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;

/**
 * Path-parallel execution with a controllable thread count. Work items write only to their own path index,
 * random numbers are derived from the path index, and all aggregations are done afterwards in path order,
 * so results do not depend on the thread count or scheduling.
 */
public final class ParallelPaths {

    private ParallelPaths() {
    }

    public static void forEach(int paths, int threads, IntConsumer action) {
        if (threads <= 1) {
            for (int p = 0; p < paths; p++) {
                action.accept(p);
            }
            return;
        }
        ForkJoinPool pool = new ForkJoinPool(threads);
        try {
            pool.submit(() -> IntStream.range(0, paths).parallel().forEach(action)).join();
        } finally {
            pool.shutdown();
        }
    }
}

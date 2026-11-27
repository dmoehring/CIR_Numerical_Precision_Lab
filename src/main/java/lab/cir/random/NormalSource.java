package lab.cir.random;

import java.util.SplittableRandom;

/**
 * Deterministic source of standard normal variates, shared by both numerical variants (common random numbers).
 *
 * <p>Generator: every path {@code p} has its own {@link SplittableRandom} (SplitMix64, stable algorithm since
 * Java 8) seeded with {@code mix64(seed + GOLDEN * (p + 1))}. A uniform {@code u = (k + 0.5) 2^-52} with the top
 * 52 bits {@code k} of {@code nextLong()} is exactly representable and lies strictly inside (0,1) and is mapped to a normal by
 * {@link InverseNormal} (AS 241). Path {@code p} therefore always receives the same normals, independently
 * of the number of paths, the evaluation order or thread scheduling; nothing has to be stored.
 *
 * <p>Brownian coupling: the source produces {@code normalsPerStep} fine normals per simulation step.
 * A coarse Brownian increment is {@code sqrt(dt / m) * (Z_1 + ... + Z_m)}; grids with
 * {@code stepsPerYear * normalsPerStep = const} thus share the same underlying fine Brownian path.
 */
public final class NormalSource {

    private static final long GOLDEN = 0x9E3779B97F4A7C15L;
    private static final double TWO_POW_MINUS_52 = 0x1.0p-52;

    private final long seed;
    private final int normalsPerPath;

    public NormalSource(long seed, int normalsPerPath) {
        if (normalsPerPath <= 0) {
            throw new IllegalArgumentException("normalsPerPath must be positive");
        }
        this.seed = seed;
        this.normalsPerPath = normalsPerPath;
    }

    public long seed() { return seed; }

    public int normalsPerPath() { return normalsPerPath; }

    /** Fills {@code out[0..normalsPerPath)} with the normals of path {@code pathIndex}. */
    public void fillPath(int pathIndex, double[] out) {
        if (out.length < normalsPerPath) {
            throw new IllegalArgumentException("buffer too small");
        }
        SplittableRandom rng = new SplittableRandom(pathSeed(seed, pathIndex));
        for (int k = 0; k < normalsPerPath; k++) {
            out[k] = InverseNormal.quantile(uniform(rng.nextLong()));
        }
    }

    public double[] path(int pathIndex) {
        double[] z = new double[normalsPerPath];
        fillPath(pathIndex, z);
        return z;
    }

    static double uniform(long bits) {
        return ((bits >>> 12) + 0.5) * TWO_POW_MINUS_52; // exact, in [2^-53, 1 - 2^-53]
    }

    static long pathSeed(long seed, int pathIndex) {
        return mix64(seed + GOLDEN * (pathIndex + 1L));
    }

    /** SplitMix64 finaliser (Steele, Lea, Flood 2014). */
    static long mix64(long z) {
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }
}

package lab.cir.model;

/**
 * Handling of the Feller condition {@code 2 * theta * mu >= sigma^2}.
 *
 * <p>The Feller condition is a property of the continuous CIR process (the origin is not attained).
 * It does <em>not</em> prevent negative values of the Euler discretisation; those are handled
 * separately by the full-truncation scheme and counted during simulation.
 */
public enum FellerMode {
    /** Default: parameters violating the Feller condition are rejected at construction time. */
    STRICT,
    /** Explicit opt-in for robustness experiments: violations are allowed (but still reported). */
    ALLOW_VIOLATION
}

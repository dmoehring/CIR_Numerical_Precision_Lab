package lab.cir.induction;

/** Summation of cluster targets in the double variant (hybrid control experiment). */
public enum SummationMode {
    NAIVE,
    /** Kahan compensated summation. */
    KAHAN
}

package lab.cir.induction;

/** Which training paths enter the continuation estimate at an exercise date. */
public enum PathSelection {
    /** All training paths (default). Exercise itself still requires a positive immediate exercise value. */
    ALL,
    /** Only paths with {@code E_j > 0} (classic Longstaff-Schwartz choice). */
    IN_THE_MONEY
}

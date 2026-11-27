package lab.cir.model;

/** Thrown in {@link FellerMode#STRICT} mode when {@code 2 * theta * mu < sigma^2}. */
public final class FellerConditionViolationException extends IllegalArgumentException {
    private static final long serialVersionUID = 1L;

    public FellerConditionViolationException(String message) {
        super(message);
    }
}

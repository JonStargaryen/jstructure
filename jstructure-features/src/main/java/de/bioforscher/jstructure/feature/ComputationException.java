package de.bioforscher.jstructure.feature;

/**
 * The exception case for feature providers.
 * Created by bittrich on 5/17/17.
 */
public class ComputationException extends RuntimeException {
    public ComputationException() {
        super();
    }

    public ComputationException(String s) {
        super(s);
    }

    public ComputationException(String s, Throwable throwable) {
        super(s, throwable);
    }

    public ComputationException(Throwable throwable) {
        super(throwable);
    }
}

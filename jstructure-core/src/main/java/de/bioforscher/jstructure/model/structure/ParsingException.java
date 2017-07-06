package de.bioforscher.jstructure.model.structure;

/**
 * Occurs when {@link ProteinParser} struggles.
 * Created by bittrich on 1/19/17.
 */
public class ParsingException extends RuntimeException {
    public ParsingException() {
        super();
    }

    public ParsingException(String s) {
        super(s);
    }

    public ParsingException(String s, Throwable throwable) {
        super(s, throwable);
    }

    public ParsingException(Throwable throwable) {
        super(throwable);
    }
}
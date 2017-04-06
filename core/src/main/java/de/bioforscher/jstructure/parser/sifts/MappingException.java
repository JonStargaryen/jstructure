package de.bioforscher.jstructure.parser.sifts;

/**
 * Occurs when to SIFTS-mapping fails for given residue numbers.
 * Created by bittrich on 4/6/17.
 */
public class MappingException extends RuntimeException {
    public MappingException() {
        super();
    }

    public MappingException(String message) {
        super(message);
    }

    public MappingException(String message, Throwable cause) {
        super(message, cause);
    }

    public MappingException(Throwable cause) {
        super(cause);
    }
}

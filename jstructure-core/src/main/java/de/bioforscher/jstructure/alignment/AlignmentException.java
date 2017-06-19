package de.bioforscher.jstructure.alignment;

/**
 * Occurs for failed alignments, e.g. when no matching atom pairs could be determined.
 * Created by bittrich on 6/19/17.
 */
public class AlignmentException extends RuntimeException {
    public AlignmentException() {
    }

    public AlignmentException(String message) {
        super(message);
    }

    public AlignmentException(String message, Throwable cause) {
        super(message, cause);
    }

    public AlignmentException(Throwable cause) {
        super(cause);
    }

    public AlignmentException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }
}

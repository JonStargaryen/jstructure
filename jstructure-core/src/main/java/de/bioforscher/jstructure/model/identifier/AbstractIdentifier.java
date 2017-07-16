package de.bioforscher.jstructure.model.identifier;

/**
 * The abstract implementation of {@link Identifier}.
 * Created by S on 14.07.2017.
 */
public abstract class AbstractIdentifier implements Identifier {
    /**
     * Enforce a meaningful toString() implementation.
     */
    public abstract String toString();

    /**
     * Enforce a meaningful equals() implementation.
     */
    public abstract boolean equals(Object other);

    /**
     * Enforce a meaningful hashCode() implementation.
     */
    public abstract int hashCode();
}

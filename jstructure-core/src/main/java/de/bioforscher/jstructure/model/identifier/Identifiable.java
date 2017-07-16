package de.bioforscher.jstructure.model.identifier;

/**
 * Specifies the capability of a container to return a specific name.
 * Created by S on 18.11.2016.
 */
public interface Identifiable<I extends Identifier> {
    /**
     * A more-or-less unique identifier.
     * @return the identifier reported by this entity
     */
    I getIdentifier();
}

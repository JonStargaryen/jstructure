package de.bioforscher.jstructure.model;

/**
 * Defines capabilities to create a copy of itself. The class must define a copy constructor.
 * Created by S on 14.07.2017.
 */
public interface Copyable {
    /**
     * Create a deep copy of this instance. {@link de.bioforscher.jstructure.model.feature.FeatureContainer} will be
     * ignored by contract.
     * @return a copy with all references (to children and parent) also cloned and set
     */
    Copyable createDeepCopy();

    /**
     * Create a shallow copy of this instance. {@link de.bioforscher.jstructure.model.feature.FeatureContainer} will be
     * ignored by contract.
     * @return a copy with all fields cloned and children and parent reference pointing to either a empty list or null
     */
    Copyable createShallowCopy();
}

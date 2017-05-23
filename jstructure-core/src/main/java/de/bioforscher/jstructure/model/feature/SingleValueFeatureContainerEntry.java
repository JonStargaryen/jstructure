package de.bioforscher.jstructure.model.feature;

/**
 * Most entries wrap a single value. This is their interface.
 * Created by bittrich on 5/23/17.
 */
public interface SingleValueFeatureContainerEntry<T> {
    T getValue();
}

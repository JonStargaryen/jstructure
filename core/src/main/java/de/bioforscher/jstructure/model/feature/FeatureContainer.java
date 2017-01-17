package de.bioforscher.jstructure.model.feature;

import java.util.List;
import java.util.Map;

/**
 * Specifies the capabilities to store arbitrary information linked to an implementing container. Map entries are
 * identified by enums.
 * Created by S on 02.10.2016.
 */
public interface FeatureContainer {
    /**
     * Handle to the internal feature map of the container. Access should be realized using the high-level functions.
     * @return the map storing all additional data attached to this element
     * @see FeatureContainer#getFeature(Class, String)
     * @see FeatureContainer#setFeature(String, Object)
     */
    Map<String, Object> getFeatureMap();

    /**
     * Retrieves an arbitrary object from the container. The expected content type has to be provided, so the returned
     * element can be casted.
     * @param expectedContentType the class of the retrieved element
     * @param key the unique identifier of the map entry
     * @param <ContentType> the class the entry will be cast to
     * @return the retrieved, casted entry
     */
    default <ContentType> ContentType getFeature(Class<ContentType> expectedContentType, String key) {
        return expectedContentType.cast(getFeatureMap().get(key));
    }

    /**
     * Retrieves an arbitrary object from the container. The expected content type has to be provided, so the returned
     * element can be casted.
     * @param expectedListContentType the class of the retrieved collection
     * @param key the unique identifier of the map entry
     * @param <ContentType> the class the entries of this list will be cast to
     * @return a collection of retrieved entries
     * @see FeatureContainer#getFeature(Class, String)
     */
    @SuppressWarnings("unchecked")
    default <ContentType> List<ContentType> getFeatureAsList(Class<ContentType> expectedListContentType, String key) {
        List entry = (List) getFeatureMap().get(key);
        if(!entry.isEmpty()) {
            // some fail-fast safety net, as the correct cast completely depends on the way how and when this method is
            // invoked, so ClassCastExceptions could arise later on
            expectedListContentType.cast(entry.get(0));
        }
        return (List<ContentType>) entry;
    }

    /**
     * Convenience function to retrieve int features.
     * @param key the unique identifier of the entry
     * @return the value as primitive int
     * @see FeatureContainer#getFeature(Class, String)
     */
    default int getFeatureAsInt(String key) {
        return getFeature(Integer.class, key);
    }

    /**
     * Convenience function to retrieve double features.
     * @param key the unique identifier of the entry
     * @return the value as primitive double
     * @see FeatureContainer#getFeature(Class, String)
     */
    default double getFeatureAsDouble(String key) {
        return getFeature(Double.class, key);
    }

    /**
     * Stores one object in the feature container.
     * @param key the enum representing this feature, old entries will be replaced
     * @param featureValue the actual content
     */
    default void setFeature(String key, Object featureValue) {
        getFeatureMap().put(key, featureValue);
    }
}

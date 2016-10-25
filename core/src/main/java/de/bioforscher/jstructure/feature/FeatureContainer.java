package de.bioforscher.jstructure.feature;

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
     * @see FeatureContainer#getFeature(Class, Enum)
     * @see FeatureContainer#setFeature(Enum, Object)
     */
    Map<Enum, Object> getFeatureMap();

    /**
     * Retrieves an arbitrary object from the container. The expected content type has to be provided, so the returned
     * element can be casted.
     * @param expectedContentType the class of the retrieved element
     * @param enumKey the unique identifier of the map entry
     * @param <ContentType> the class the entry will be cast to
     * @return the retrieved, casted entry
     */
    default <ContentType> ContentType getFeature(Class<ContentType> expectedContentType, Enum enumKey) {
        return expectedContentType.cast(getFeatureMap().get(enumKey));
    }

    /**
     * Convenience function to retrieve int features.
     * @param enumKey the unique identifier of the entry
     * @return the value as primitive int
     * @see FeatureContainer#getFeature(Class, Enum)
     */
    default int getIntFeature(Enum enumKey) {
        return getFeature(Integer.class, enumKey);
    }

    /**
     * Convenience function to retrieve double features.
     * @param enumKey the unique identifier of the entry
     * @return the value as primitive double
     * @see FeatureContainer#getFeature(Class, Enum)
     */
    default double getDoubleFeature(Enum enumKey) {
        return getFeature(Double.class, enumKey);
    }

    /**
     * Stores one object in the feature container.
     * @param enumKey the enum representing this feature, old entries will be replaced
     * @param featureValue the actual content
     */
    default void setFeature(Enum enumKey, Object featureValue) {
        getFeatureMap().put(enumKey, featureValue);
    }
}

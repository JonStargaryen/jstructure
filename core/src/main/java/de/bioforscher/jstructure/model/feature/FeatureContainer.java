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
     * Retrieves an arbitrary object from the container. The expected content type has to be provided, so the returned
     * element can be casted.
     * @param expectedListContentType the class of the retrieved collection
     * @param enumKey the unique identifier of the map entry
     * @param <ContentType> the class the entries of this list will be cast to
     * @return a collection of retrieved entries
     * @see FeatureContainer#getFeature(Class, Enum)
     */
    @SuppressWarnings("unchecked")
    default <ContentType> List<ContentType> getFeatureAsList(Class<ContentType> expectedListContentType, Enum enumKey) {
        List entry = (List) getFeatureMap().get(enumKey);
        if(!entry.isEmpty()) {
            //some fail-fast safety net, as the correct cast completely depends on the way how this method is invoked
            // and ClassCastException could arise later on
            expectedListContentType.cast(entry.get(0));
        }
        return (List<ContentType>) entry;
    }

    /**
     * Convenience function to retrieve int features.
     * @param enumKey the unique identifier of the entry
     * @return the value as primitive int
     * @see FeatureContainer#getFeature(Class, Enum)
     */
    default int getFeatureAsInt(Enum enumKey) {
        return getFeature(Integer.class, enumKey);
    }

    /**
     * Convenience function to retrieve double features.
     * @param enumKey the unique identifier of the entry
     * @return the value as primitive double
     * @see FeatureContainer#getFeature(Class, Enum)
     */
    default double getFeatureAsDouble(Enum enumKey) {
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

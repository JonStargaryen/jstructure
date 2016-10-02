package de.bioforscher.jstructure.feature;

import java.util.Map;

/**
 * Specifies the capabilities to store arbitrary information linked to an implementing container.
 * TODO: is Map<String, Object> the best possible way? enums would forbid any external extensions - Strings are somewhat unsafe though - maybe {@link FeatureProvider} should enumerate the feature names they are assigning
 * Created by S on 02.10.2016.
 */
public interface FeatureContainer {
    /**
     *
     * @return
     */
    Map<String, Object> getFeatureMap();

    /**
     *
     * @param contentType
     * @param featureKey
     * @param <ContentType>
     * @return
     */
    default <ContentType> ContentType getFeature(Class<ContentType> contentType, String featureKey) {
        return contentType.cast(getFeatureMap().get(featureKey));
    }

    /**
     *
     * @param featureKey
     * @return
     */
    default double getDoubleFeature(String featureKey) {
        return getFeature(double.class, featureKey);
    }

    /**
     *
     * @param featureKey
     * @param value
     */
    default void setFeature(String featureKey, Object value) {
        getFeatureMap().put(featureKey, value);
    }
}

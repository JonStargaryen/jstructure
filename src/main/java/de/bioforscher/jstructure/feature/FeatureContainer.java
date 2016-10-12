package de.bioforscher.jstructure.feature;

import java.util.Map;

/**
 * Specifies the capabilities to store arbitrary information linked to an implementing container.
 * Created by S on 02.10.2016.
 */
public interface FeatureContainer {
    /**
     *
     * @return
     */
    Map<String, Object> getFeatureMap();

    /**
     * Retrieves one entry from the feature map.
     * @param featureKey the name of the element to retrieve
     * @param <ContentType> the content of this feature
     * @return the value of the requested feature or <code>null</code>
     */
    default <ContentType> ContentType getFeature(Class<ContentType> expectedContentType, String featureKey) {
        return expectedContentType.cast(getFeatureMap().get(featureKey));
    }

    /**
     *
     * @param expectedContentType
     * @param enumKey
     * @param <ContentType>
     * @return
     */
    default <ContentType> ContentType getFeature(Class<ContentType> expectedContentType, Enum enumKey) {
        return getFeature(expectedContentType, enumKey.name());
    }

    default double getDoubleFeature(String featureKey) {
        return getFeature(double.class, featureKey);
    }

    default double getFeature(Enum enumKey) {
        return getDoubleFeature(enumKey.name());
    }

    default int getIntFeature(String featureKey) {
        return getFeature(int.class, featureKey);
    }

    default int getIntFeature(Enum enumKey) {
        return getIntFeature(enumKey.name());
    }

    /**
     * Stores one object in the feature container.
     * @param featureKey the name of this feature, old entries will be replaced
     * @param featureValue the actual content
     */
    default void setFeature(String featureKey, Object featureValue) {
        getFeatureMap().put(featureKey, featureValue);
    }

    /**
     *
     * @param enumKey
     * @param featureValue
     */
    default void setFeature(Enum enumKey, Object featureValue) {
        setFeature(enumKey.name(), featureValue);
    }
}

package de.bioforscher.jstructure.model.feature;

import java.util.HashSet;
import java.util.Set;

/**
 * The root instance for entities who itself or whose children are {@link Featureable} instances. Provide
 * capabilities to determine whether certain features are present or not.
 * Created by bittrich on 5/17/17.
 */
@Deprecated
public interface FeatureContainerRoot extends Featureable {
    default boolean featureIsPresent(Class<? extends FeatureContainerEntry> contentClass) {
        //TODO need tests
        return getFeatureContainer()
                .getFeatureOptional(GlobalFeatureContainer.class)
                .get()
                .featureIsPresent(contentClass);
    }

    default void registerFeature(Class<? extends FeatureContainerEntry> contentClass) {
        getFeatureContainer()
                .getFeatureOptional(GlobalFeatureContainer.class)
                .get()
                .registerFeature(contentClass);
    }

    @Deprecated
    class GlobalFeatureContainer extends FeatureContainerEntry {
        private Set<Class<? extends FeatureContainerEntry>> registeredEntries;

        GlobalFeatureContainer() {
            // no feature provider is associated to this synthetic instance
            super(null);
            this.registeredEntries = new HashSet<>();
        }

        boolean featureIsPresent(Class<? extends FeatureContainerEntry> contentClass) {
            return registeredEntries.stream()
                    .anyMatch(contentClass::isAssignableFrom);
        }

        void registerFeature(Class<? extends FeatureContainerEntry> contentClass) {
            registeredEntries.add(contentClass);
        }
    }
}

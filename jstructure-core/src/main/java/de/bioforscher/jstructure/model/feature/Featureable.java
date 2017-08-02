package de.bioforscher.jstructure.model.feature;

/**
 * An entity providing access to {@link FeatureContainer}.
 * Created by S on 28.04.2017.
 */
public interface Featureable {
    FeatureContainer getFeatureContainer();

    <C extends FeatureContainerEntry> C getFeature(Class<C> contentClass);

    void setFeatureContainer(FeatureContainer featureContainer);
}

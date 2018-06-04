package de.bioforscher.jstructure.model.feature;

/**
 * The abstract implementation of an entity capable of providing access to a feature container. Basically, the only
 * reasoning behind this class is to keep the implementation of other classes rather clean and separate the
 * functionality of the container from other methods directly linked to an entity.
 * Created by bittrich on 5/17/17.
 */
public abstract class AbstractFeatureable implements Featureable {
    private FeatureContainer featureContainer;

    protected AbstractFeatureable() {
        this.featureContainer = new FeatureContainer();
    }

    @Override
    public FeatureContainer getFeatureContainer() {
        return featureContainer;
    }

    @Override
    public void setFeatureContainer(FeatureContainer featureContainer) {
        this.featureContainer = featureContainer;
    }
}

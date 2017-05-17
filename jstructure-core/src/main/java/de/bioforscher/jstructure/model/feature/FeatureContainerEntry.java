package de.bioforscher.jstructure.model.feature;

/**
 * An entry of the {@link FeatureContainer}. Provides specific information of the entry and meta-information of the
 * computation (respectively the {@link de.bioforscher.jstructure.model.feature.FeatureProvider} which was used to
 * compute it. By contract stored values are expected to be final.
 * Created by S on 28.04.2017.
 */
public class FeatureContainerEntry {
    private final AbstractFeatureProvider featureProvider;

    public FeatureContainerEntry(AbstractFeatureProvider featureProvider) {
        this.featureProvider = featureProvider;
    }

    /**
     * Access to the feature provider which was used to compute this entry.
     * @return the instance which was used to compute this entry
     */
    public AbstractFeatureProvider getFeatureProvider() {
        return featureProvider;
    }
}

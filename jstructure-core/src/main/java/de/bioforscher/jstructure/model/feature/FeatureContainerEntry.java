package de.bioforscher.jstructure.model.feature;

/**
 * An entry of the {@link FeatureContainer}. Provides specific information of the entry and meta-information of the
 * computation (respectively the {@link de.bioforscher.jstructure.model.feature.FeatureProvider} which was used to
 * calculate it. By contract stored values are expected to be final.
 * Created by S on 28.04.2017.
 */
public class FeatureContainerEntry {
    private final AbstractFeatureProvider featureProvider;

    /**
     * The standard constructor for new feature map entries.
     * @param featureProvider the provider which created this feature - may be <code>null</code>, however it is strongly
     *                        advised to provide this information for full-fledged {@link AbstractFeatureProvider}
     *                        implementations (the use case it have different implementations annotating the same
     *                        feature by different means, by providing this generator information entries of the same
     *                        content class can be distinguished)
     */
    public FeatureContainerEntry(AbstractFeatureProvider featureProvider) {
        this.featureProvider = featureProvider;
    }

    /**
     * Access to the feature provider which was used to calculate this entry.
     * @return the instance which was used to calculate this entry
     */
    public AbstractFeatureProvider getFeatureProvider() {
        return featureProvider;
    }
}

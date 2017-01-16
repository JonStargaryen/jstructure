package de.bioforscher.jstructure.model.feature;

/**
 * The abstract implementation of each {@link FeatureProvider}. Implicitly registers each provider in a global
 * "registry" implemented by {@link FeatureServiceRegistry}.
 * Created by S on 16.01.2017.
 */
public abstract class AbstractFeatureProvider implements FeatureProvider {
    public AbstractFeatureProvider() {
        registerFeatureProvider();
    }

    private void registerFeatureProvider() {
        FeatureServiceRegistry.getInstance().registerFeatureProvider(this);
    }
}

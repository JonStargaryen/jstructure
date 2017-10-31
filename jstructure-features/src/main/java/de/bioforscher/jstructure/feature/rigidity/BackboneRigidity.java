package de.bioforscher.jstructure.feature.rigidity;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

@DefaultFeatureProvider(DynaMineBridge.class)
public class BackboneRigidity extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double backboneRigidity;

    public BackboneRigidity(FeatureProvider featureProvider, double backboneRigidity) {
        super(featureProvider);
        this.backboneRigidity = backboneRigidity;
    }

    public double getBackboneRigidity() {
        return backboneRigidity;
    }

    @Override
    public Double getValue() {
        return getBackboneRigidity();
    }
}

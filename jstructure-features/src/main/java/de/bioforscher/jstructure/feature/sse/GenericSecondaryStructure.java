package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * The wrapper for secondary structure elements.
 * Created by bittrich on 6/28/17.
 */
public class GenericSecondaryStructure extends FeatureContainerEntry {
    protected SecondaryStructureElement secondaryStructure;

    public GenericSecondaryStructure(AbstractFeatureProvider featureProvider, SecondaryStructureElement secondaryStructure) {
        super(featureProvider);
        this.secondaryStructure = secondaryStructure;
    }

    public SecondaryStructureElement getSecondaryStructure() {
        return secondaryStructure;
    }

    public void setSecondaryStructure(SecondaryStructureElement secondaryStructure) {
        this.secondaryStructure = secondaryStructure;
    }
}

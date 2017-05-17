package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

import java.util.List;

/**
 * Represents the membrane layer of proteins as a set of 2 layers of points.
 * Created by bittrich on 5/17/17.
 */
public class GenericMembrane extends FeatureContainerEntry {
    private List<double[]> membraneAtoms;

    public GenericMembrane(AbstractFeatureProvider featureProvider, List<double[]> membraneAtoms) {
        super(featureProvider);
        this.membraneAtoms = membraneAtoms;
    }

    public List<double[]> getMembraneAtoms() {
        return membraneAtoms;
    }
}

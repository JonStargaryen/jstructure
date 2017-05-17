package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

import java.util.List;

/**
 * The collection of trans-membrane helices.
 * Created by bittrich on 5/17/17.
 */
public class TransMembraneHelixContainer extends FeatureContainerEntry {
    private List<GenericTransMembraneHelix> transMembraneHelices;

    public TransMembraneHelixContainer(AbstractFeatureProvider featureProvider) {
        super(featureProvider);

    }

    public List<GenericTransMembraneHelix> getTransMembraneHelices() {
        return transMembraneHelices;
    }

    public void addTransMembraneHelix(GenericTransMembraneHelix transMembraneHelix) {
        transMembraneHelices.add(transMembraneHelix);
    }

    //TODO checks if some residue is trans-membrane or not
}

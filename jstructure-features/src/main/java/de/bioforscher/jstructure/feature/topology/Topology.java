package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * Indicates the topology of this amino acid.
 * Created by bittrich on 6/8/17.
 */
@DefaultFeatureProvider(OrientationsOfProteinsInMembranesAnnotator.class)
public class Topology extends FeatureContainerEntry {
    private final boolean transmembrane;

    public Topology(FeatureProvider featureProvider, boolean transmembrane) {
        super(featureProvider);
        this.transmembrane = transmembrane;
    }

    public boolean isTransmembrane() {
        return transmembrane;
    }
}

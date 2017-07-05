package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * Indicates the topology of this amino acid.
 * Created by bittrich on 6/8/17.
 */
public class Topology extends FeatureContainerEntry {
    private final boolean transmembrane;

    public Topology(AbstractFeatureProvider featureProvider, boolean transmembrane) {
        super(featureProvider);
        this.transmembrane = transmembrane;
    }

    public boolean isTransmembrane() {
        return transmembrane;
    }
}

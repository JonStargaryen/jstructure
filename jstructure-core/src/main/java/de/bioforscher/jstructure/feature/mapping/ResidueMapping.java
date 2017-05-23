package de.bioforscher.jstructure.feature.mapping;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * Residue-specific mappings, e.g. to the UniProt sequence number of particular residues.
 * Created by bittrich on 5/17/17.
 */
public class ResidueMapping extends FeatureContainerEntry {
    public static final String UNKNOWN_MAPPING = "?";
    private final String uniProtResidueNumber;

    ResidueMapping(AbstractFeatureProvider featureProvider, String uniProtResidueNumber) {
        super(featureProvider);
        this.uniProtResidueNumber = uniProtResidueNumber;
    }

    public ResidueMapping(AbstractFeatureProvider featureProvider) {
        this(featureProvider, UNKNOWN_MAPPING);
    }

    public String getUniProtResidueNumber() {
        return uniProtResidueNumber;
    }
}

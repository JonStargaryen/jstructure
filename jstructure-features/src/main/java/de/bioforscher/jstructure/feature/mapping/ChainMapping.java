package de.bioforscher.jstructure.feature.mapping;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * Chain-specific mappings, e.g. to EC numbers for chains.
 * Created by bittrich on 5/17/17.
 */
public class ChainMapping extends FeatureContainerEntry {
    public static final String UNKNOWN_MAPPING = "?";
    private final String uniProtId;
    private final String ecNumber;
    private final String pfamId;

    ChainMapping(AbstractFeatureProvider featureProvider, String uniProtId, String ecNumber, String pfamId) {
        super(featureProvider);
        this.uniProtId = uniProtId;
        this.ecNumber = ecNumber;
        this.pfamId = pfamId;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public String getEcNumber() {
        return ecNumber;
    }

    public String getPfamId() {
        return pfamId;
    }
}

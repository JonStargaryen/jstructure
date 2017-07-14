package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.Map;

/**
 * Extract a sequence conversation profile from a multiple-sequence alignment.
 * Created by bittrich on 7/13/17.
 */
public interface SequenceConservationCalculator {
    /**
     * Extract a sequence conservation profile for a given chain
     * @param alignment an alignment contain this chain's sequence
     * @param chain the chain to process
     */
    void extractConservationProfile(Map<String, String> alignment, Chain chain);
}

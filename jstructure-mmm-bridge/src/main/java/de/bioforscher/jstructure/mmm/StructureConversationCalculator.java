package de.bioforscher.jstructure.mmm;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.mmm.model.Itemset;

import java.util.List;
import java.util.Map;

/**
 * Compose a profile describing the degree of structural conservation of all amino acids.
 * Created by bittrich on 7/13/17.
 */
public interface StructureConversationCalculator {
    /**
     * Annotate a protein with a structural conservation profile.
     * @param extractedItemsets mmm results
     * @param chain the chain of that the profile should be computed of
     */
    List<Double> extractConservationProfile(Map<Itemset<String>, List<Itemset<String>>> extractedItemsets, Chain chain);
}

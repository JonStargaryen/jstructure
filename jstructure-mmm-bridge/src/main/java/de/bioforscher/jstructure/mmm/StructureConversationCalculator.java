package de.bioforscher.jstructure.mmm;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.mmm.model.Itemset;

import java.util.List;
import java.util.Map;

/**
 * Compose a profile describing the degree of structural conservation of all amino acids.
 * Created by bittrich on 7/13/17.
 */
public interface StructureConversationCalculator {
    /**
     * Annotate a protein with {@link StructureConservationProfile} entries.
     * @param extractedItemsets mmm results
     * @param protein the protein the profile could be assigned to
     */
    void extractConservationProfile(Map<Itemset<String>, List<Itemset<String>>> extractedItemsets, Structure protein);
}

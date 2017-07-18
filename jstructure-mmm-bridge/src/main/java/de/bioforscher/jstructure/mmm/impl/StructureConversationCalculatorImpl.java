package de.bioforscher.jstructure.mmm.impl;

import de.bioforscher.jstructure.mmm.StructureConversationCalculator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.mmm.model.Itemset;
import de.bioforscher.singa.chemistry.physical.model.LeafIdentifier;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * The implemenation of the structural conservation calculator.
 * Created by bittrich on 7/13/17.
 */
public class StructureConversationCalculatorImpl implements StructureConversationCalculator {
    @Override
    public List<Double> extractConservationProfile(Map<Itemset<String>, List<Itemset<String>>> extractedItemsets, Chain chain) {
        return chain.aminoAcids()
                .map(aminoAcid -> computeConservationScore(extractedItemsets, aminoAcid))
                .collect(Collectors.toList());
    }

    private double computeConservationScore(Map<Itemset<String>, List<Itemset<String>>> extractedItemsets, AminoAcid aminoAcid) {
        String pdbId = aminoAcid.getParentChain().getParentStructure().getProteinIdentifier().getPdbId();
        String chainId = aminoAcid.getParentChain().getChainIdentifier().getChainId();
        int residueNumber = aminoAcid.getResidueIdentifier().getResidueNumber();

        return extractedItemsets.values()
                .stream()
                .flatMap(Collection::stream)
                .map(Itemset::getStructuralMotif)
                .filter(Optional::isPresent)
                .map(Optional::get)
                // check if this itemset describes the given amino acid
                .filter(structuralMotif -> structuralMotif.getLeafSubstructures()
                        .stream()
                        .anyMatch(leafSubstructure -> {
                            LeafIdentifier identifier = leafSubstructure.getLeafIdentifier();
                            return identifier.getPdbIdentifier().equalsIgnoreCase(pdbId) &&
                                    identifier.getChainIdentifer().equals(chainId) &&
                                    identifier.getIdentifier() == residueNumber;
                        }))
                .count() / (double) extractedItemsets.size();
    }
}

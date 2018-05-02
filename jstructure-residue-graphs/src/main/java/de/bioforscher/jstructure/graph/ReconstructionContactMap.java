package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureType;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a protein as contact map which can be used for reconstruction. Ignores covalent bonds and local contacts.
 */
public class ReconstructionContactMap extends ResidueGraph {
    private final String sequence;
    private final String secondaryStructureElements;
    private String name;

    public static ReconstructionContactMap createReconstructionContactMap(Chain chain) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        return new ReconstructionContactMap(aminoAcids, createContactList(aminoAcids, InteractionScheme.CALPHA8));
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public ReconstructionContactMap(List<AminoAcid> aminoAcids, List<Pair<AminoAcid, AminoAcid>> contacts) {
        super(aminoAcids, filterForLongRangeContacts(contacts));
        this.sequence = aminoAcids.stream()
                .map(AminoAcid::getOneLetterCode)
                .collect(Collectors.joining());
        this.secondaryStructureElements = aminoAcids.stream()
                .map(aminoAcid -> aminoAcid.getFeature(GenericSecondaryStructure.class))
                .map(GenericSecondaryStructure::getSecondaryStructure)
                .map(SecondaryStructureType::getReducedRepresentation)
                .map(String::toUpperCase)
                .collect(Collectors.joining());
    }

    private static List<Pair<AminoAcid, AminoAcid>> filterForLongRangeContacts(List<Pair<AminoAcid, AminoAcid>> contacts) {
        return contacts.stream()
                .filter(ResidueGraph::isLongRangeContact)
                .collect(Collectors.toList());
    }

    public String getSequence() {
        return ">CHAIN" + System.lineSeparator() + sequence;
    }

    public String getSecondaryStructureElements() {
        return secondaryStructureElements;
    }

    public String getCaspRRRepresentation() {
        // i  j  d1  d2  p
        // 24 33 0 8 12.279515
        // i and j: residue numbers
        return contacts.stream()
                .map(edge -> (edge.getLeft().getAminoAcidIndex() + 1) + " " +
                        //TODO changeable contact definition
                        (edge.getRight().getAminoAcidIndex() + 1) + " 0 8.0 1")
                .collect(Collectors.joining(System.lineSeparator(),
                        sequence + System.lineSeparator(),
                        ""));
    }

    public int getNumberOfAminoAcids() {
        return aminoAcids.size();
    }

    public int getNumberOfContacts() {
        return contacts.size();
    }
}

package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The container object of a found sequence motif.
 * Created by S on 02.10.2016.
 */
public class SequenceMotif {
    private final SequenceMotifDefinition motifDefinition;
    private final ChainIdentifier chainId;
    private final int startResidueNumber;
    private final int endResidueNumber;
    private final List<AminoAcid> aminoAcids;

    SequenceMotif(SequenceMotifDefinition candidate,
                  ChainIdentifier chainId,
                  int startResidueNumber,
                  int endResidueNumber,
                  List<AminoAcid> aminoAcids) {
        this.motifDefinition = candidate;
        this.chainId = chainId;
        this.startResidueNumber = startResidueNumber;
        this.endResidueNumber = endResidueNumber;
        this.aminoAcids = aminoAcids;
    }

    public SequenceMotifDefinition getMotifDefinition() {
        return motifDefinition;
    }

    public ChainIdentifier getChainId() {
        return chainId;
    }

    public int getStartResidueNumber() {
        return startResidueNumber;
    }

    public int getEndResidueNumber() {
        return endResidueNumber;
    }

    public String getSequence() {
        return aminoAcids.stream()
                .map(AminoAcid::getOneLetterCode)
                .collect(Collectors.joining());
    }

    public List<AminoAcid> getAminoAcids() {
        return aminoAcids;
    }

    @Override
    public String toString() {
        return motifDefinition + " " + chainId + "_" + startResidueNumber + "-" + endResidueNumber;
    }
}
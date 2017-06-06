package de.bioforscher.jstructure.feature.motif;

/**
 * The container object of a found sequence motif.
 * Created by S on 02.10.2016.
 */
public class SequenceMotif {
    private final SequenceMotifDefinition motifDefinition;
    private final String chainId;
    private final int startResidueNumber;
    private final int endResidueNumber;
    private final String sequence;

    SequenceMotif(SequenceMotifDefinition candidate, String chainId, int startResidueNumber, int endResidueNumber, String sequence) {
        this.motifDefinition = candidate;
        this.chainId = chainId;
        this.startResidueNumber = startResidueNumber;
        this.endResidueNumber = endResidueNumber;
        this.sequence = sequence;
    }

    public SequenceMotifDefinition getMotifDefinition() {
        return motifDefinition;
    }

    public String getChainId() {
        return chainId;
    }

    public int getStartResidueNumber() {
        return startResidueNumber;
    }

    public int getEndResidueNumber() {
        return endResidueNumber;
    }

    public String getSequence() {
        return sequence;
    }
}
package de.bioforscher.jstructure.feature.motif;

/**
 * The container object of a found sequence motif.
 * Created by S on 02.10.2016.
 */
public class SequenceMotif {
    private final SequenceMotifDefinition motifDefinition;
    private final String chainId;
    private final int start;
    private final int end;
    private final String sequence;

    public SequenceMotif(SequenceMotifDefinition candidate, String chainId, int start, int end, String sequence) {
        this.motifDefinition = candidate;
        this.chainId = chainId;
        this.start = start;
        this.end = end;
        this.sequence = sequence;
    }

    public SequenceMotifDefinition getMotifDefinition() {
        return motifDefinition;
    }

    public String getChainId() {
        return chainId;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String getSequence() {
        return sequence;
    }
}
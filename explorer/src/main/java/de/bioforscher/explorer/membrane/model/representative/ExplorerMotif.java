package de.bioforscher.explorer.membrane.model.representative;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;

/**
 * The reduced representation of {@link SequenceMotif} objects.
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerMotif {
    private String chain, type, sequence;
    private int start, end;

    public ExplorerMotif() {

    }

    public ExplorerMotif(SequenceMotif sequenceMotif) {
        this.chain = sequenceMotif.getStartGroup().getParentChain().getChainId();
        this.type = sequenceMotif.getMotifDefinition().name();
        this.start = sequenceMotif.getStartGroup().getResidueNumber();
        this.end = sequenceMotif.getEndGroup().getResidueNumber();
        this.sequence = sequenceMotif.getGroupContainer().getAminoAcidSequence();
    }

    public String getChain() {
        return chain;
    }

    public String getType() {
        return type;
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

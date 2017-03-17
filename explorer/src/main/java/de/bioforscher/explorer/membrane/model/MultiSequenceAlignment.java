package de.bioforscher.explorer.membrane.model;

/**
 * Represent the multi-sequence alignemnt of 1 cluster.
 * Created by bittrich on 3/17/17.
 */
public class MultiSequenceAlignment {
    private String representativeChainId, alignment;

    public MultiSequenceAlignment() {
    }

    public MultiSequenceAlignment(String representativeChainId, String alignment) {
        this.representativeChainId = representativeChainId;
        this.alignment = alignment;
    }

    public String getRepresentativeChainId() {
        return representativeChainId;
    }

    public String getAlignment() {
        return alignment;
    }
}

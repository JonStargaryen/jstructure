package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;

import java.util.Map;
import java.util.Set;

/**
 * The implementation of a multiple-sequence alignment result.
 * Created by S on 12.07.2017.
 */
class MultipleSequenceAlignmentResultImpl implements MultipleSequenceAlignmentResult {
    private final Map<String, String> alignment;

    MultipleSequenceAlignmentResultImpl(Map<String, String> alignment) {
        this.alignment = alignment;
    }

    @Override
    public Set<String> getIdsOfAlignedSequences() {
        return alignment.keySet();
    }

    @Override
    public Map<String, String> getAlignedSequences() {
        return alignment;
    }
}

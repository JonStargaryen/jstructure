package de.bioforscher.jstructure.align;

import java.util.List;

/**
 * Specification of a service providing multiple sequence alignments.
 * Created by bittrich on 7/12/17.
 */
public interface MultipleSequenceAligner {
    /**
     * Computes a multiple-sequence alignment for a given collection of sequences.
     * @param fastaSequences the input sequences - expected to be protein sequence in FASTA format
     * @return the results wrapped as {@link MultipleSequenceAlignmentResult}
     * @throws AlignmentException when no alignment could be achieved
     */
    MultipleSequenceAlignmentResult align(List<String> fastaSequences) throws AlignmentException;
}

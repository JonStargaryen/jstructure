package de.bioforscher.jstructure.align;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

/**
 * Specification of a service providing multiple sequence alignments.
 * Created by bittrich on 7/12/17.
 */
public interface MultipleSequenceAligner {
    /**
     * Computes a multiple-sequence alignment for a given collection of sequences.
     * @param fastaSequences the input sequences - expected to be protein sequence in FASTA format
     * @return a map with the id as key and the aligned sequence as value
     * @throws ExecutionException upon failed execution
     */
    Map<String, String> process(List<String> fastaSequences) throws ExecutionException;
}

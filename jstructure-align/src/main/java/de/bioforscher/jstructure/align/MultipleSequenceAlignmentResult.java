package de.bioforscher.jstructure.align;

import java.util.Map;
import java.util.Set;

/**
 * The result of a multiple-sequence alignment.
 * Created by S on 12.07.2017.
 */
public interface MultipleSequenceAlignmentResult {
    /**
     * Returns the identifiers of all aligned sequences as collection of strings.
     * @return a collection of all ids
     */
    Set<String> getIdsOfAlignedSequences();

    /**
     * The alignment as map with ids as keys and the aligned sequence as value. All sequence have equal length and may
     * contain gaps denoted as '-'.
     * @return a map (id -> alignedSequence)
     */
    Map<String, String> getAlignedSequences();
}

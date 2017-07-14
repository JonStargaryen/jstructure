package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;

import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Shared functions by both clustal omega 'implementations'.
 * Created by bittrich on 7/14/17.
 */
abstract class AbstractClustalOmega implements MultipleSequenceAligner {
    private static final Pattern SEQUENCE_PATTERN = Pattern.compile(">");
    private static final Pattern LINE_PATTERN = Pattern.compile("\\s+");

    MultipleSequenceAlignmentResult createMultipleSequenceAlignmentResult(String rawAlignmentString) {
        Map<String, String> alignment = new HashMap<>();
        SEQUENCE_PATTERN.splitAsStream(rawAlignmentString)
                // skip first (empty) pair
                .skip(1)
                .forEach(line -> {
                    // split alignment at newlines
                    String[] split = LINE_PATTERN.split(line);
                    // substring to drop FASTA-character
                    alignment.put(split[0],
                            // skip id and join all other lines
                            Stream.of(split).skip(1).collect(Collectors.joining()));
                });

        return new MultipleSequenceAlignmentResultImpl(alignment);
    }
}

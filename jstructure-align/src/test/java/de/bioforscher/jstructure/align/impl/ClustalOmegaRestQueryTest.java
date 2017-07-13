package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Tests for a MSA by ClustalOmega.
 * Created by S on 12.07.2017.
 */
public class ClustalOmegaRestQueryTest {
    private MultipleSequenceAligner multipleSequenceAligner;
    private List<String> sequences;

    @Before
    public void setup() {
        multipleSequenceAligner = new ClustalOmegaRestQuery();
        Protein protein1 = ProteinParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        Protein protein2 = ProteinParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1BRR))
                .minimalParsing(true)
                .parse();
        sequences = Stream.of(protein1, protein2)
                .flatMap(Protein::chainsWithAminoAcids)
                .map(chain -> ">" + chain.getChainIdentifier() + System.lineSeparator() + chain.getAminoAcidSequence())
                .collect(Collectors.toList());
    }

    @Test
    public void shouldExecuteClustalOmegaQuery() {
        MultipleSequenceAlignmentResult alignment = multipleSequenceAligner.align(sequences);
        // assert all sequences are equally long
        Assert.assertEquals("aligned sequences do not match in length", 1, alignment.getAlignedSequences()
                .values()
                .stream()
                .map(String::length)
                .distinct()
                .count());
    }
}

package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
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
        Structure protein1 = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        Structure protein2 = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1BRR))
                .minimalParsing(true)
                .parse();
        sequences = Stream.of(protein1, protein2)
                .flatMap(Structure::chainsWithAminoAcids)
                .map(chain -> ">" + chain.getChainIdentifier() + System.lineSeparator() + chain.getAminoAcidSequence())
                .collect(Collectors.toList());
    }

    @Test
    @Ignore
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

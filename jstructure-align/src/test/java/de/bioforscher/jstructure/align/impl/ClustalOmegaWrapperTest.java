package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Test for the wrapper of clustal omega.
 * Created by bittrich on 7/14/17.
 */
public class ClustalOmegaWrapperTest {
    private List<String> sequences;
    private ClustalOmegaWrapper clustalOmegaWrapper;

    @Before
    public void setup() {
        Random random = new Random();
        sequences = IntStream.range(0, 10)
                .mapToObj(i -> ">sequence_" + (i + 1) + System.lineSeparator() +
                    IntStream.range(0, 20 + random.nextInt(20))
                        .mapToObj(p -> AminoAcid.Family.values()[random.nextInt(20)].getOneLetterCode())
                        .collect(Collectors.joining()))
                .collect(Collectors.toList());
        this.clustalOmegaWrapper = new ClustalOmegaWrapper();
    }

    @Test
    public void mapShouldFeatureReferenceAsFirstEntry() {
        Map<String, String> alignmentMap = clustalOmegaWrapper.align(sequences).getAlignedSequences();
        Assert.assertEquals(sequences.size(), alignmentMap.size());
        Assert.assertEquals("first entry is not reference",
                "sequence_1",
                alignmentMap.entrySet().stream()
                .findFirst()
                .get()
                .getKey());
    }

    @Test
    public void shouldAlignSequences() {
        Map<String, String> alignmentMap = clustalOmegaWrapper.align(sequences).getAlignedSequences();
        Assert.assertEquals(sequences.size(), alignmentMap.size());
        alignmentMap.entrySet().forEach(System.out::println);
    }
}
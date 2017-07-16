package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Test for the sequence conservation annotator.
 * Created by bittrich on 7/12/17.
 */
@Deprecated
public class SequenceConservationAnnotatorTest {
    private Chain chain;
    private Map<String, String> alignmentMap;

    @Before
    public void setup() {
        chain = MutationJobImpl.createProtein("Q", "ABCDEFGHIKLMNPQRSTVWY").getChains().get(0);
        alignmentMap = new HashMap<>();
        alignmentMap.put("Q", "ABCDEFGHIKLMNPQRSTVWY");
        alignmentMap.put("1", "AB-DEFGHIKLMNPQRSTVWY");
        alignmentMap.put("2", "AB-DFFGHIKLMNPQRSTV--");
        alignmentMap.put("3", "AB-DFFGHIKLMMPQRSAV--");
        alignmentMap.put("4", "ABCDEFGHIKL-NPQRSTVW-");
        alignmentMap.put("5", "GBCDEFGHIKL-MPQRSTVWY");
        alignmentMap.put("6", "AB-DGFGHIKL-NPQRSAVW-");
        alignmentMap.put("7", "GB-DEFGHIKL-NPQRSTVW-");
        alignmentMap.put("8", "ABCDEFGHIKL-MPQRSTVWY");
        alignmentMap.put("9", "ABCDGFGHIKL-NPQRSAVWY");
    }

    @Test
    public void shouldAnnotateConservationProfile() {
        SequenceConservationAnnotator.process(chain, alignmentMap);
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        aminoAcids.forEach(aminoAcid -> {
            Optional<SequenceConservationProfile> sequenceConservationProfile = aminoAcid.getFeatureContainer().getFeatureOptional(SequenceConservationProfile.class);
            Assert.assertTrue("sequence conservation profile not present for " + aminoAcid, sequenceConservationProfile.isPresent());
            System.out.println(sequenceConservationProfile.get());
        });
    }
}
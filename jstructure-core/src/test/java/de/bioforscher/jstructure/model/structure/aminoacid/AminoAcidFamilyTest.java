package de.bioforscher.jstructure.model.structure.aminoacid;

import org.junit.Assert;
import org.junit.Test;

/**
 * Test for the Family enum of AminoAcid.
 * Created by bittrich on 7/12/17.
 */
public class AminoAcidFamilyTest {
    @Test
    public void shouldResolveOneLetterCodeCorrectly() {
        // olc 'M' is used by methionine and selenomethionine
        Assert.assertEquals("resolving of one-letter-codes corrupted",
                AminoAcid.Family.METHIONINE,
                AminoAcid.Family.resolveOneLetterCode("M"));
    }
}
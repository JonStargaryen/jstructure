package de.bioforscher.jstructure.model.structure.identifier;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit tests for ProteinIdentifier.
 * Created by bittrich on 4/27/17.
 */
public class PdbIdTest {
    @Test(expected = IllegalArgumentException.class)
    public void shouldFailForInvalidPdbId() {
        ProteinIdentifier.createFromPdbId("xxxx");
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailForAdditionalNameContainingUnderscore() {
        ProteinIdentifier.createFromAdditionalName("1bbr_test");
    }

    @Test
    public void shouldConvertIdToLowerCase() {
        Assert.assertEquals("1brr", ProteinIdentifier.createFromPdbId("1BRR").getPdbId());
    }

    @Test
    public void shouldCreateComplexIdSeparatedByHyphen() {
        Assert.assertTrue(ProteinIdentifier.createFromPdbIdAndName("1brr", "mutated").getFullName().contains("-"));
    }
}
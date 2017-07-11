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
        IdentifierFactory.createProteinIdentifier("xxxx");
    }

    @Test
    public void shouldConvertIdToLowerCase() {
        Assert.assertEquals("1brr", IdentifierFactory.createProteinIdentifier("1BRR").getPdbId());
    }

    @Test
    public void shouldCreateComplexIdSeparatedByHyphen() {
        Assert.assertTrue(IdentifierFactory.createProteinIdentifier("1brr", "mutated").getFullName().contains("-"));
    }
}
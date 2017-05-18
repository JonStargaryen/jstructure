package de.bioforscher.jstructure.model.structure.identifier;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit tests for PdbId.
 * Created by bittrich on 4/27/17.
 */
public class PdbIdTest {
    @Test(expected = IllegalArgumentException.class)
    public void shouldFailForInvalidPdbId() {
        PdbId.createFromPdbId("xxxx");
    }

    @Test
    public void shouldConvertIdToLowerCase() {
        Assert.assertEquals("1brr", PdbId.createFromPdbId("1BRR").getPdbId());
    }

    @Test
    public void shouldCreateComplexIdSeparatedByHyphen() {
        Assert.assertTrue(PdbId.createFromPdbIdAndName("1brr", "mutated").getFullName().contains("-"));
    }
}
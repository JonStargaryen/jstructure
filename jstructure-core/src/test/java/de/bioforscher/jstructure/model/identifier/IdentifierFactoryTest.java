package de.bioforscher.jstructure.model.identifier;

import org.junit.Assert;
import org.junit.Test;

/**
 * Test for the identifier factory and its validation of ids.
 * Created by S on 14.07.2017.
 */
public class IdentifierFactoryTest {
    @Test
    public void shouldResolveOldUniProtId() {
        Assert.assertEquals("UniProt id was not resolved correctly",
                "P38554",
                IdentifierFactory.createUniProtIdentifier("CYC32_DESNO").getUniProtId());
    }

    @Test
    public void shouldCreateChainId() {
        Assert.assertEquals("ChainIdentifier format corrupted",
                "1brr_A",
                IdentifierFactory.createChainIdentifier("1brr", "A").getFullName());
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailForInvalidPdbId() {
        IdentifierFactory.createProteinIdentifier("xxxx");
    }

    @Test
    public void shouldConvertIdToLowerCase() {
        Assert.assertEquals("ProteinIdentifier format corrupted",
                "1brr",
                IdentifierFactory.createProteinIdentifier("1BRR").getPdbId());
    }

    @Test
    public void shouldCreateComplexIdSeparatedByHyphen() {
        Assert.assertTrue("ComplexName format corrupted",
                IdentifierFactory.createProteinIdentifier("1brr",
                "mutated").getFullName().contains("-"));
    }
}
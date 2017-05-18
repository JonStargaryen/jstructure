package de.bioforscher.jstructure.model.structure.identifier;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit tests for PdbChainId.
 * Created by bittrich on 5/18/17.
 */
public class PdbChainIdTest {
    @Test
    public void shouldCreateChainId() {
        Assert.assertEquals("1brr_A", PdbChainId.createFromChainId(PdbId.createFromPdbId("1brr"), "A").getFullName());
    }
}
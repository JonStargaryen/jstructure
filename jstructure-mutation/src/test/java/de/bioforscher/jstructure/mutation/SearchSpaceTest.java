package de.bioforscher.jstructure.mutation;

import org.junit.Assert;
import org.junit.Test;

/**
 * Test for the search space provider.
 * Created by bittrich on 7/6/17.
 */
public class SearchSpaceTest {
    @Test
    public void shouldGetIdsOfWholePdb() {
        Assert.assertTrue("full pdb search space seems corrupted", SearchSpace.fullPdb().size() > 130000);
    }

    @Test
    public void shouldGetIdsOfPF00959() {
        // lysozyme
        SearchSpace.pdbByPfam("PF00959").forEach(System.out::println);
    }

    @Test
    public void shouldGetIdsOfP00720() {
        // lysozyme
        SearchSpace.pdbByUniProt("P00720").forEach(System.out::println);
    }
}
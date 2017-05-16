package de.bioforscher.jstructure.parser.pdb;

import org.junit.Assert;
import org.junit.Test;

import java.util.List;

/**
 * Test the PDB REST interface.
 * Created by bittrich on 3/17/17.
 */
public class PDBDatabaseQueryTest {
    @Test
    public void shouldFetchClusterFor1m0l() {
        List<String> cluster = PDBDatabaseQuery.fetchSequenceCluster("1m0l", "A");
        Assert.assertTrue(cluster.size() >= 118);
        Assert.assertTrue(cluster.contains("1BRR.A"));
        Assert.assertTrue(cluster.contains("1BRR.B"));
        Assert.assertTrue(cluster.contains("1BRR.C"));
        Assert.assertTrue(cluster.contains("1MOL.A"));
    }
}
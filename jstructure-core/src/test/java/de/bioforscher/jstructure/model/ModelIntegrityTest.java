package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Checks especially the createCopy() methods of the data model.
 * Created by S on 23.11.2016.
 */
public class ModelIntegrityTest {
    private Protein protein;

    @Before
    public void setup() {
        this.protein = ProteinParser.source("1brr").parse();
    }

    @Test
    public void testGroupCreateCopy() {
        AtomContainer groupOriginal = protein.getGroups().get(0);
        AtomContainer groupCopy = groupOriginal.createCopy();
        Assert.assertTrue(groupCopy instanceof Group);
    }

    @Test
    public void testChainCreateCopy() {
        GroupContainer chainCopy = protein.getChains()
                .get(0)
                .createCopy();
        Assert.assertTrue(chainCopy instanceof Chain);
    }

    @Test
    public void testProteinCreateCopy() {
        ChainContainer proteinCopy = protein.createCopy();
        Assert.assertTrue(proteinCopy instanceof Protein);
    }
}

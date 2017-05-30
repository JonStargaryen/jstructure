package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test capabilities of the structure collectors.
 * Created by bittrich on 5/29/17.
 */
public class StructureCollectorsTest {
    private Protein protein;

    @Before
    public void setup() {
        protein = ProteinParser.source("1brr").parse();
    }

    @Test
    public void toAtomContainer() throws Exception {
        //TODO rework concept, has there to be some structure for atoms (i.e. organize selection again into chains e.g.)
        AtomContainer container = protein.atoms()
                .collect(StructureCollectors.toAtomContainer());
        container.atoms().forEach(atom -> {
            Assert.assertEquals("parent reference was lost",
                    "1brr",
                    atom.getParentGroup().getParentChain().getParentProtein().getPdbId().getPdbId());
        });
    }

    @Test
    public void toGroupContainer() throws Exception {
    }

    @Test
    public void toChainContainer() throws Exception {
    }

}
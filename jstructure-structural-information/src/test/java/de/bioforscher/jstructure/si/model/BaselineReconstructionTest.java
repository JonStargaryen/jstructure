package de.bioforscher.jstructure.si.model;

import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinition;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class BaselineReconstructionTest {
    private Chain chain;
    private ContactDefinition contactDefinition;

    @Before
    public void setup() {
        chain = StructureParser.fromInputStream(TestUtils.getResourceAsInputStream("confold/short.pdb"))
                .parse()
                .getFirstChain();
        contactDefinition = ContactDefinitionFactory.createAlphaCarbonContactDefinition(8.0);
    }

    @Test
    public void testComputeQ() {
        ReconstructionContactMap reference = ReconstructionContactMap.createReconstructionContactMap(chain, contactDefinition);
        ReconstructionContactMap reconstruct = ReconstructionContactMap.createReconstructionContactMap(
                StructureParser.fromInputStream(TestUtils.getResourceAsInputStream("confold/short.pdb"))
                        .parse()
                        .getFirstChain(),
                contactDefinition);
        Assert.assertEquals("for the same contact map Q should be 1",
                1.0,
                BaselineReconstruction.computeQ(reference, reconstruct),
                TestUtils.TOLERANT_ERROR_MARGIN);
    }
}
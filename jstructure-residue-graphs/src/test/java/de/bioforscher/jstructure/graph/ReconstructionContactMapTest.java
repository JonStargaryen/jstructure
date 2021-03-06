package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Test;

public class ReconstructionContactMapTest {
    @Test
    public void shouldReportCorrectSecondaryStructureElementString() {
        Chain chain = StructureParser.fromPdbId("1mbc")
                .parse()
                .getFirstChain();
        System.out.println(ReconstructionContactMap.createReconstructionContactMap(chain,
                ContactDefinitionFactory.createAlphaCarbonContactDefinition(8.0)).getSecondaryStructureElements());
    }
}
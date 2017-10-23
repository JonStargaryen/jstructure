package de.bioforscher.jstructure.feature.rigidity;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Test;

public class DynaMineBridgeTest {
    @Test
    public void shouldPredictBackboneRigidity() {
        Structure structure = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        new DynaMineBridge().process(structure);
    }
}
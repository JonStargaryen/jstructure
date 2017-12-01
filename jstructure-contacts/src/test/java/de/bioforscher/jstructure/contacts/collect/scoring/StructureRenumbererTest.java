package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Test;

public class StructureRenumbererTest {
    @Test
    public void shouldRenumberSequences() {
        Structure reference = StructureParser.source(TestUtils.getResourceAsInputStream("casp/T0515.pdb")).parse();
        Structure model = StructureParser.source(TestUtils.getResourceAsInputStream("casp/T0515TS001_3")).parse();
        StructureRenumberer.renumberStructure(reference, model);
    }
}
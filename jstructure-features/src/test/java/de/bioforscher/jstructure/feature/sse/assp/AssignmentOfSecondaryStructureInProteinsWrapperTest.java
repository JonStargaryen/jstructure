package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Test;

/**
 * Tests the CLI-wrapper for ASSP.
 * Created by bittrich on 7/5/17.
 */
public class AssignmentOfSecondaryStructureInProteinsWrapperTest {
    @Test
    public void shouldExecuteASSP() {
        Structure protein = StructureParser.fromPdbId("2jho").parse();
        new AssignmentOfSecondaryStructureInProteinsWrapper().process(protein);
        protein.aminoAcids()
                .forEach(aminoAcid -> System.out.println(aminoAcid.getIdentifier() + " " +
                        aminoAcid.getFeature(GenericSecondaryStructure.class)
                                .getSecondaryStructure()
                                .getOneLetterRepresentation()));
    }
}
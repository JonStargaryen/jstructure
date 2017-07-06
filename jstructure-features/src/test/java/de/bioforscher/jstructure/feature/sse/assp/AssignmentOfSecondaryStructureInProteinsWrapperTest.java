package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Test;

/**
 * Tests the CLI-wrapper for ASSP.
 * Created by bittrich on 7/5/17.
 */
public class AssignmentOfSecondaryStructureInProteinsWrapperTest {
    @Test
    public void shouldExecuteASSP() {
        Protein protein = ProteinParser.source("2jho").parse();
        new AssignmentOfSecondaryStructureInProteinsWrapper().process(protein);
        protein.aminoAcids()
                .forEach(aminoAcid -> System.out.println(aminoAcid.getIdentifier() + " " +
                        aminoAcid.getFeatureContainer()
                                .getFeature(GenericSecondaryStructure.class)
                                .getSecondaryStructure()
                                .getOneLetterRepresentation()));
    }
}
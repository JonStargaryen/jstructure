package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

/**
 * Functional test for the ASSP implementation.
 * Created by bittrich on 6/28/17.
 * @see <a href="http://nucleix.mbu.iisc.ac.in/assp/outputfiles.html">http://nucleix.mbu.iisc.ac.in/assp/outputfiles.html</a>
 *      for additional information on the nature of output files
 */
public class AssignmentOfSecondaryStructureInProteinsTest {
    private static final AssignmentOfSecondaryStructureInProteins assignmentOfSecondaryStructureInProteins = new AssignmentOfSecondaryStructureInProteins();

    @Test
    public void shouldAgreeWithReferenceOutput() {
        Protein protein = ProteinParser.source("2jho")
                .minimalParsing(true)
                .parse();

        assignmentOfSecondaryStructureInProteins.process(protein);

        //TODO test for agreement with 2jho_nh.out, 2jho_cont.out and 2jho_assp.out
        protein.aminoAcids()
                .forEach(aminoAcid -> {
                    ASSPSecondaryStructure asspSecondaryStructure = aminoAcid.getFeatureContainer().getFeature(ASSPSecondaryStructure.class);
                    System.out.println(asspSecondaryStructure.getAlpha() + " " +
                            asspSecondaryStructure.getThree() + " " +
                            asspSecondaryStructure.getPi() + " " +
                            String.valueOf(aminoAcid.getResidueNumber().getResidueNumber()));
                });
    }
}
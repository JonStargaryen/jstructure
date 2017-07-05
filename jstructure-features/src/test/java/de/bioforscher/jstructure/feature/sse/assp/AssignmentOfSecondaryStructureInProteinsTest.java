package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Functional test for the ASSP implementation.
 * Created by bittrich on 6/28/17.
 * @see <a href="http://nucleix.mbu.iisc.ac.in/assp/outputfiles.html">http://nucleix.mbu.iisc.ac.in/assp/outputfiles.html</a>
 *      for additional information on the nature of output files
 */
public class AssignmentOfSecondaryStructureInProteinsTest {
    private static final AssignmentOfSecondaryStructureInProteins assignmentOfSecondaryStructureInProteins = new AssignmentOfSecondaryStructureInProteins();
    private Protein protein;
    private List<String[]> contLines;
    private List<String[]> asspLines;

    @Before
    public void setup() throws IOException {
        protein = ProteinParser.source("2jho")
                .minimalParsing(true)
                .parse();

        assignmentOfSecondaryStructureInProteins.process(protein);

        contLines = TestUtils.getResourceAsStream("assp/2jho_cont.out")
                .map(string -> string.replaceAll("\\s+", " "))
                .map(string -> string.split(" "))
                .filter(split -> split.length == 18)
                .collect(Collectors.toList());

        asspLines = TestUtils.getResourceAsStream("assp/2jho_assp.out")
                .map(string -> string.replaceAll("\\s+", " "))
                .map(string -> string.split(" "))
                .collect(Collectors.toList());
    }

    @Test
    public void shouldAgreeInAssignment() {
        // ASSP assigns 3 states to each amino acid in a stretch
        protein.aminoAcids()
                .forEach(aminoAcid -> {
                    String chainId = aminoAcid.getParentChain().getChainId().getChainId();
                    int residueNumber = aminoAcid.getResidueNumber().getResidueNumber();
                    ASSPSecondaryStructure asspSecondaryStructure = aminoAcid.getFeatureContainer().getFeature(ASSPSecondaryStructure.class);

                    contLines.stream()
                            .filter(split -> split[12].equals(chainId))
                            .filter(split -> Integer.valueOf(split[4]) == residueNumber)
                            .forEach(split -> {
                                if(!split[0].equals(asspSecondaryStructure.getAlpha()) ||
                                        !split[1].equals(asspSecondaryStructure.getThree()) ||
                                        !split[2].equals(asspSecondaryStructure.getPi())) {
//                                    System.err.println("assignment does not match for " + aminoAcid.getIdentifier() +
//                                            " expected: " + split[0] + " " + split[1] + " " + split[2] + " found: " +
//                                            asspSecondaryStructure);
                                    Assert.fail("assignment does not match for " + aminoAcid.getIdentifier() +
                                            " expected: " + split[0] + " " + split[1] + " " + split[2] + " found: " +
                                            asspSecondaryStructure);
                                }
                            });
                });

        //TODO deviations for
        // ALPHA_HELIX	21	32	12    vs    AlphaHelix   2 VAL A   21  SER A   35  1                                  15
        // ALPHA_HELIX	33	36	4
    }
}
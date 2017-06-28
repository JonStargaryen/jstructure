package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucState;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Checks whether the ported BioJava DSSP implementation is still in agreement with the original one.
 * Created by S on 01.11.2016.
 */
public class DictionaryOfProteinSecondaryStructureTest {
    private static final String ID = "1brr";
    private DictionaryOfProteinSecondaryStructure featureProvider;

    @Before
    public void setup() {
        featureProvider = new DictionaryOfProteinSecondaryStructure();
    }

    @Test
    public void shouldClearPseudoAtomsAfterRun() {
        Protein protein = ProteinParser.source(ID).parse();
        featureProvider.process(protein);
        boolean containsPseudoHydrogenLine = protein.getPdbRepresentation().contains("ATOM      0  H");
        Assert.assertFalse("pseudo-atoms were not removed!", containsPseudoHydrogenLine);
    }

    @Test
    public void testTorsionAngleComputationForResiduesInDifferentChains() {
        // should ignore amino acids in different chains
        featureProvider.process(ProteinParser.source("4cqn").parse());
    }

    @Test
    public void checkAgreement() throws IOException, StructureException {
        checkAgreement(ID);
    }

    private void checkAgreement(String id) {
        try {
            String jstructureAnnotation = getSecondaryStructureAnnotation(id);
            String biojavaAnnotation = getDSSPAnnotatedStructure(id);

            Assert.assertEquals(biojavaAnnotation, jstructureAnnotation);
        } catch (Exception e) {
            Assert.fail(e.getMessage());
        }
    }

    @Test
    public void checkAgreementForErrorCases() throws IOException, StructureException {
        Stream.of("1ar1", "3aqp", "4b4a")
                .forEach(this::checkAgreement);
    }

    private String getSecondaryStructureAnnotation(String id) {
        // load structure
        Protein protein = ProteinParser.source(id).parse();
        // assign states
        featureProvider.process(protein);

        // return complete DSSP annotation string from jstructrue
        return protein.aminoAcids()
                .map(residue -> residue.getFeatureContainer().getFeature(DSSPSecondaryStructure.class))
                .map(DSSPSecondaryStructure::getSecondaryStructure)
                .map(SecondaryStructureElement::getOneLetterRepresentation)
                .map(character -> character.equals("c") ? " " : character)
                .collect(Collectors.joining());
    }

    private String getDSSPAnnotatedStructure(String id) throws IOException, StructureException {
        // load structure
        Structure protein = new PDBFileReader().getStructureById(id);
        // assign states
        new SecStrucCalc().calculate(protein, true);

        // return complete DSSP annotation string from BioJava
        return protein.getChains()
                      .stream()
                      .flatMap(chain -> chain.getAtomGroups(GroupType.AMINOACID).stream())
                      .map(aminoAcid -> aminoAcid.getProperty(Group.SEC_STRUC))
                      .map(SecStrucState.class::cast)
                      .map(SecStrucState::getType)
                      .map(type -> String.valueOf(type.type))
                      .collect(Collectors.joining());
    }
}

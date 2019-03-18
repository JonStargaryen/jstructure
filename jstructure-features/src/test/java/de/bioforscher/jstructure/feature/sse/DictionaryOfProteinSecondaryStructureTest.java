package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
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
        Structure protein = StructureParser.fromPdbId(ID).parse();
        featureProvider.process(protein);
        boolean containsPseudoHydrogenLine = protein.getPdbRepresentation().contains("ATOM      0  H");
        Assert.assertFalse("pseudo-atoms were not removed!", containsPseudoHydrogenLine);
    }

    @Test
    public void testTorsionAngleComputationForResiduesInDifferentChains() {
        // should ignore amino acids in different chains
        featureProvider.process(StructureParser.fromPdbId("4cqn").parse());
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
        Structure protein = StructureParser.fromPdbId(id).parse();
        // assign states
        featureProvider.process(protein);

        // return complete DSSP annotation string from jstructrue
        return protein.aminoAcids()
                .map(residue -> residue.getFeature(DSSPSecondaryStructure.class))
                .map(DSSPSecondaryStructure::getSecondaryStructure)
                .map(SecondaryStructureType::getOneLetterRepresentation)
                .map(character -> character.equals("c") ? " " : character)
                .collect(Collectors.joining());
    }

    private String getDSSPAnnotatedStructure(String id) throws IOException, StructureException {
        // load structure
        org.biojava.nbio.structure.Structure protein = new PDBFileReader().getStructureById(id);
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

    @Test
    public void test1bta() {
        Structure structure = StructureParser.fromPdbId("1bta").parse();
        featureProvider.process(structure);
        structure.aminoAcids()
                .forEach(aminoAcid -> {
            System.out.println(aminoAcid);
            DSSPSecondaryStructure sse = aminoAcid.getFeature(DSSPSecondaryStructure.class);
            System.out.println(sse.getSecondaryStructure().name());
        });

        System.out.println(structure.getFirstChain()
                .aminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(DSSPSecondaryStructure.class).getSecondaryStructure().getOneLetterRepresentation())
                //.map(c -> !c.equals(" ") ? c : "-")
                .collect(Collectors.joining()));
    }

    @Test
    public void test1btaBiojava() throws IOException, StructureException {
        org.biojava.nbio.structure.Structure structure = new PDBFileReader().getStructureById("1bta");
        new SecStrucCalc().calculate(structure, true);

        // return complete DSSP annotation string from BioJava
        System.out.println(structure.getChains()
                .stream()
                .flatMap(chain -> chain.getAtomGroups(GroupType.AMINOACID).stream())
                .map(aminoAcid -> aminoAcid.getProperty(Group.SEC_STRUC))
                .map(SecStrucState.class::cast)
                .map(SecStrucState::getType)
                .map(type -> String.valueOf(type.type))
                .map(type -> type.equals(" ") ? "-" : type)
                .collect(Collectors.joining()));
    }
}

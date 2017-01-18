package feature.sse;

import de.bioforscher.jstructure.feature.sse.DSSPSecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
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

/**
 * Checks whether the ported BioJava DSSP implementation is still in agreement with the original one.
 * Created by S on 01.11.2016.
 */
public class SecondaryStructureAnnotatorFunctionalTest {
    private static final String ID = "1brr";
    private AbstractFeatureProvider featureProvider;

    @Before
    public void setup() {
        featureProvider = FeatureProviderRegistry.getInstance().resolve(SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES);
    }

    @Test
    public void checkAgreement() throws IOException, StructureException {
        String id = ID;

        String jstructureAnnotation = getSecondaryStructureAnnotation(id);
        String biojavaAnnotation = getDSSPAnnotatedStructure(id);

        //TODO all sheets/bridges are skipped by jstructure - otherwise these are in agreement even though the test validly fails at the moment
        Assert.assertEquals(biojavaAnnotation, jstructureAnnotation);
    }

    private String getSecondaryStructureAnnotation(String id) {
        // load structure
        Protein protein = ProteinParser.parseProteinById(ID);
        // assign states
        featureProvider.process(protein);

        // return complete DSSP annotation string from jstructrue
        return Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getFeature(de.bioforscher.jstructure.feature.sse.SecStrucState.class,
                        SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES))
                .map(de.bioforscher.jstructure.feature.sse.SecStrucState::getSecondaryStructure)
                .map(DSSPSecondaryStructureElement::getOneLetterRepresentation)
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

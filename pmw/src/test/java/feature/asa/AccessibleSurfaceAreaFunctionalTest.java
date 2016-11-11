package feature.asa;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.structure.asa.GroupAsa;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Checks for agreement between BioJava's ASA calculator and the jstructure adaptation.
 * Created by S on 03.11.2016.
 */
public class AccessibleSurfaceAreaFunctionalTest {
    private static final String ID = "1acj";

    @Test
    public void checkAgreement() throws IOException, StructureException {
        String id = ID;

        List<Double> jstructureASA = getJStructureASA(id);
        List<Double> biojavaASA = getBioJavaASA(id);

        //TODO implement real test, respectively fix differences
        Pair.sequentialPairsOf(jstructureASA, biojavaASA).forEach(doublePair ->
            Assert.assertEquals(doublePair.getLeft(), doublePair.getRight(), 0.001)
        );
    }

    private static List<Double> getJStructureASA(String id) {
        // load structure
        Protein protein = ProteinParser.parseProteinById(id);
        // assign states
        new AccessibleSurfaceAreaCalculator().process(protein);

        // return complete DSSP annotation string from jstructrue
        return protein.residues()
                .map(residue ->
                        residue.getDoubleFeature(AccessibleSurfaceAreaCalculator.FeatureNames.ACCESSIBLE_SURFACE_AREA))
                .collect(Collectors.toList());
    }

    private static List<Double> getBioJavaASA(String id) throws IOException, StructureException {
        // load structure
        Structure protein = new PDBFileReader().getStructureById(id);

        // assign ASA
        return Arrays.stream(new AsaCalculator(protein,
                     AsaCalculator.DEFAULT_PROBE_SIZE,
                     AsaCalculator.DEFAULT_N_SPHERE_POINTS,
                     AsaCalculator.DEFAULT_NTHREADS,
                     false).getGroupAsas())
                     .map(GroupAsa::getAsaU)
                     .collect(Collectors.toList());
    }
}

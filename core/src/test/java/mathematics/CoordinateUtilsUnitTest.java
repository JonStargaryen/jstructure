package mathematics;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;

import static util.TestUtils.TOLERANT_ERROR_MARGIN;

/**
 * Unit tests for CoordinateUtils.
 * Created by S on 24.10.2016.
 */
public class CoordinateUtilsUnitTest {
    private Protein protein1acj;
    private Protein protein1brr;

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParser.parseProteinById("1acj");
        protein1brr = ProteinParser.parseProteinById("1brr");
    }

    @Test
    public void testAtomPairsInContact() {
        final double cutoff = 8.0;
        protein1acj.atomPairsInContact(cutoff)
                   .filter(pair -> LinearAlgebra3D.distance(pair.getFirst().getCoordinates(),
                                                            pair.getSecond().getCoordinates()) > cutoff)
                   .forEach(pair -> Assert.fail("pairs should be closer than " + cutoff +
                               ", this is violated for pair: " + pair));
    }

    @Test
    public void testResiduePairsInContact() {
        final double cutoff = 8.0;
        protein1acj.residuePairsInContact(cutoff)
                   .filter(pair -> LinearAlgebra3D.distance(pair.getFirst().getAlphaCarbon().getCoordinates(),
                        pair.getSecond().getAlphaCarbon().getCoordinates()) > cutoff)
                   .forEach(pair -> Assert.fail("pairs should be closer than " + cutoff +
                        ", this is violated for pair: " + pair));
    }

    @Test
    public void testCenterOfMassWithProtein() {
        double[] com1 = CoordinateUtils.centerOfMass(protein1acj);
        double[] com2 = centerOfMass(protein1acj);
        System.out.println("api-centerOfMass: " + Arrays.toString(com1) + " test-centerOfMass: " + Arrays.toString(com2));

        Assert.assertArrayEquals(com1, com2, 0.0);
    }

    /**
     * Old school of computing the center of mass for a collection of atoms.
     * @param protein the container
     * @return the center of mass as double[]
     */
    private double[] centerOfMass(Protein protein) {
        double[] total = new double[3];
        double mass = 0;

        for(Atom atom : protein.atoms().collect(Collectors.toList())) {
            total = LinearAlgebra3D.add(total, LinearAlgebra3D.multiply(atom.getCoordinates(), atom.getElement().getAtomicMass()));
            mass += atom.getElement().getAtomicMass();
        }

        return LinearAlgebra3D.divide(total, mass);
    }

    @Test
    public void testMaximalExtentWithProtein() {
        double maximalExtent = CoordinateUtils.maximalExtent(protein1acj);
        System.out.println(maximalExtent);
        //TODO this needs a real test
    }

    @Test
    public void shouldResultInCenteredStructure() {
        CoordinateUtils.center(protein1acj);
        Assert.assertArrayEquals(CoordinateUtils.centroid(protein1acj),
                AlignmentAlgorithm.NEUTRAL_TRANSLATION, TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldReturnNoDifferenceInRMSDForIdenticalProtein() {
        double rmsd = CoordinateUtils.calculateRMSD(protein1acj, protein1acj);
        Assert.assertEquals(0.0, rmsd, 0.0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailForRMSDCalculationOfDifferentProteins() {
        CoordinateUtils.calculateRMSD(protein1acj, protein1brr);
    }
}

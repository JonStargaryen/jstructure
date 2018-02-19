package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Tests for coordinate manipulations.
 * Created by bittrich on 5/23/17.
 */
public class AtomLinearAlgebraTest {
    private Structure protein;

    @Before
    public void setup() throws IOException {
        protein = StructureParser.fromInputStream(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
    }

    @Test
    public void center() throws Exception {
        protein.calculate().center();
        Assert.assertArrayEquals(LinearAlgebra.on(protein).centroid().getValue(),
                Transformation.NEUTRAL_TRANSLATION, TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void centerOfMass() throws Exception {
        Assert.assertArrayEquals("center of mass does not match expectation",
                new double[] { 4.7529, 65.5514, 56.9086 },
                protein.calculate().centerOfMass().getValue(),
                TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void centroid() throws Exception {
        Assert.assertArrayEquals("centroid does not match expectation",
                new double[] { 4.7625, 65.5635, 56.8967 },
                protein.calculate().centroid().getValue(),
                TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void maximalExtent() throws Exception {
        Assert.assertEquals("maximal extent does not match expectation",
                41.148,
                protein.calculate().maximalExtent(),
                TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void transform() throws Exception {
        // example from http://stackoverflow.com/questions/34050929/3d-point-rotation-algorithm
        double[] point1 = { 1, 0, 0 };
        final double[] expectedPoint = { 0, -1, 0 };
        System.out.printf("rotating %s by -90 degree yields ", Arrays.toString(point1));
        // describes a rotation by 90°
        double[][] rotation = {{ 0, -1, 0 },{ 1, 0, 0 },{ 0, 0, 1 }};
        Transformation transformation = new Transformation(rotation);
        Atom atom = Atom.builder(Element.C, point1).build();
        transformation.transform(atom);
        System.out.println(Arrays.toString(atom.getCoordinates()));
        Assert.assertArrayEquals(expectedPoint, atom.getCoordinates(), 0.0);

        double[] point2 = { 3, 0, 0 };
        final double[] expectedPoint2 = { -2, -3, 0 };
        double[] translation2 = { -2, 0, 0 };
        System.out.printf("translating %s by 2 length units toward origin and rotating it by -90 degree yields ",
                Arrays.toString(point2));
        // describes a rotation by 90°
        double[][] rotation2 = { { 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 1 } };
        Transformation transformation2 = new Transformation(translation2, rotation2);
        Atom atom2 = Atom.builder(Element.C, point2).build();
        transformation2.transform(atom2);
        System.out.println(Arrays.toString(atom2.getCoordinates()));
        Assert.assertArrayEquals(expectedPoint2, atom2.getCoordinates(), 0.0);
    }
}
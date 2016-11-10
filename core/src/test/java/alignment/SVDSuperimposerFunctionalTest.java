package alignment;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Tests the capabilities of the SVDSuperimposer.
 * Created by S on 02.10.2016.
 */
public class SVDSuperimposerFunctionalTest {
    private Protein protein1acj;
    private Protein protein1brr;
    private static final String FRAGMENT_DIR = "alignment/fragment/";

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParser.parseProteinById("1acj");
        protein1brr = ProteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldAlignFragments() throws IOException {
        Files.list(Paths.get(TestUtils.getResourceAsFilepath(FRAGMENT_DIR)))
             .forEach(System.out::println);
    }

    @Test
    public void shouldResultInPerfectAlignment() {
        SVDSuperimposer svd = new SVDSuperimposer();
        AlignmentResult result = svd.align(protein1acj, protein1acj);
        Assert.assertEquals(result.getRmsd(), 0.0, 0.001);
        Assert.assertArrayEquals(result.getTranslation(), AlignmentAlgorithm.NEUTRAL_TRANSLATION, 0.001);
    }

    @Test
    public void shouldAlignSimilarStructureFragments() {
        //TODO implement
    }

    @Test
    public void shouldResultInPerfectAlignmentForTransformedCopy() throws IOException {
        Protein protein1acjCopy = ProteinParser.parseProteinById("1acj");
        double[] translation = new double[] { 10, 20, 30 };
        CoordinateUtils.transform(protein1acjCopy,
                translation,
                AlignmentAlgorithm.NEUTRAL_ROTATION);

        SVDSuperimposer svd = new SVDSuperimposer();
        AlignmentResult result = svd.align(protein1acj, protein1acjCopy);
        Assert.assertEquals(result.getRmsd(), 0.0, 0.001);
        Assert.assertArrayEquals(result.getTranslation(), LinearAlgebra3D.multiply(translation, -1.0), 0.001);
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailAsAtomOrderIsDifferent() {
        SVDSuperimposer svd = new SVDSuperimposer();
        svd.align(protein1acj, protein1brr);
    }
}

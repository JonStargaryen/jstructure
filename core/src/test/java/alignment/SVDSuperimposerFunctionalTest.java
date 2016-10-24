package alignment;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * Created by S on 02.10.2016.
 */
public class SVDSuperimposerFunctionalTest {
    private ProteinParser proteinParser;
    private Protein protein1acj;
    private Protein protein1brr;

    @Before
    public void setup() throws IOException {
        proteinParser = new ProteinParser();
        protein1acj = proteinParser.parseProteinById("1acj");
        protein1brr = proteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldResultInPerfectAlignment() {
        SVDSuperimposer svd = new SVDSuperimposer();
        AlignmentResult result = svd.align(protein1acj, protein1acj);
        Assert.assertEquals(result.getRmsd(), 0.0, 0.001);
        Assert.assertArrayEquals(result.getTranslation(), AlignmentAlgorithm.NEUTRAL_TRANSLATION, 0.001);
    }

    @Test
    public void shouldResultInPerfectAlignmentForTransformedCopy() throws IOException {
        Protein protein1acjCopy = proteinParser.parseProteinById("1acj");
        CoordinateUtils.transform(protein1acjCopy.atoms(),
                new double[] { 10, 20, 30},
                AlignmentAlgorithm.NEUTRAL_ROTATION);

        SVDSuperimposer svd = new SVDSuperimposer();
        AlignmentResult result = svd.align(protein1acj, protein1acjCopy);
        Assert.assertEquals(result.getRmsd(), 0.0, 0.001);
        Assert.assertArrayEquals(result.getTranslation(), AlignmentAlgorithm.NEUTRAL_TRANSLATION, 0.001);
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailAsAtomOrderIsDifferent() {
        SVDSuperimposer svd = new SVDSuperimposer();
        svd.align(protein1acj, protein1brr);
    }
}

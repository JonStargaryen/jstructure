package calculation;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.stream.Stream;

/**
 * Created by S on 30.09.2016.
 */
public class CoordinateUtilsFunctionalTest {
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
    public void shouldResultInCenteredStructure() {
        Assert.assertArrayEquals(CoordinateUtils.centroid(CoordinateUtils.center(protein1acj.atoms())),
                AlignmentAlgorithm.NEUTRAL_TRANSLATION, 0.001);
    }

    @Test
    public void shouldReturnIdenticalTransformedCoordinates() {
        Atom atom = protein1acj.atoms().findFirst().get();
        final double[] coordinates = atom.getCoordinates();
        CoordinateUtils.transform(Stream.of(atom), AlignmentAlgorithm.NEUTRAL_TRANSLATION, AlignmentAlgorithm.NEUTRAL_ROTATION);
        Assert.assertArrayEquals(coordinates, atom.getCoordinates(), 0.0);
    }

    @Test
    public void shouldReturnNoDifferenceInRMSDForIdenticalProtein() {
        double rmsd = CoordinateUtils.calculateRMSD(protein1acj.atoms(), protein1acj.atoms());
        Assert.assertEquals(0.0, rmsd, 0.0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailForRMSDCalculationOfDifferentProteins() {
        CoordinateUtils.calculateRMSD(protein1acj.atoms(), protein1brr.atoms());
    }
}

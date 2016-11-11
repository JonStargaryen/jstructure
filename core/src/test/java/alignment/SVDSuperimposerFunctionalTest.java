package alignment;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.StructureCollectors;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import parser.ProteinParserFunctionalTest;
import util.TestUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

/**
 * Tests the capabilities of the SVDSuperimposer.
 * Created by S on 02.10.2016.
 */
public class SVDSuperimposerFunctionalTest {
    private Protein protein1acj;
    private Protein protein1brr;
    private static final String FRAGMENT_DIR = "alignment/fragment/";
    private GroupContainer container1;
    private GroupContainer container2;
    private GroupContainer container3;
    private GroupContainer container4;

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParserFunctionalTest.parseFilepath(ProteinParserFunctionalTest.PDB_DIRECTORY + "1acj" +
                ProteinParserFunctionalTest.PDB_EXTENSION);
        protein1brr = ProteinParserFunctionalTest.parseFilepath(ProteinParserFunctionalTest.PDB_DIRECTORY + "1brr" +
                ProteinParserFunctionalTest.PDB_EXTENSION);

        Residue his1 = new Residue("HIS", 1);
        Residue asp1 = new Residue("ASP", 2);
        Residue ser1 = new Residue("SER", 3);
        his1.addAtom(new Atom("CA", Element.C, new double[] { 6.994, 8.354, 42.405 }));
        asp1.addAtom(new Atom("CA", Element.C, new double[] { 9.429, 7.479, 48.266 }));
        ser1.addAtom(new Atom("CA", Element.C, new double[] { 5.547, 0.158, 42.050 }));
        container1 = GroupContainer.of(Stream.of(his1, asp1, ser1));

        Residue his2 = new Residue("HIS", 1);
        Residue asp2 = new Residue("ASP", 2);
        Residue ser2 = new Residue("SER", 3);
        his2.addAtom(new Atom("CA", Element.C, new double[] { 3.908, 12.066, -6.159 }));
        asp2.addAtom(new Atom("CA", Element.C, new double[] { 4.588, 6.531, -9.119 }));
        ser2.addAtom(new Atom("CA", Element.C, new double[] { 12.080, 12.645, -7.073 }));
        container2 = GroupContainer.of(Stream.of(his2, asp2, ser2));

        Residue ala1 = new Residue("HIS", 1);
        Residue his3 = new Residue("ASP", 2);
        Residue cys1 = new Residue("SER", 3);
        ala1.addAtom(new Atom("CA", Element.C, new double[] { 5.055, 74.865, 22.585 }));
        his3.addAtom(new Atom("CA", Element.C, new double[] { 7.320, 76.960, 20.325 }));
        cys1.addAtom(new Atom("CA", Element.C, new double[] { 6.021, 74.874, 17.385 }));
        container3 = GroupContainer.of(Stream.of(ala1, his3, cys1));

        Residue ala2 = new Residue("HIS", 1);
        Residue his4 = new Residue("ASP", 2);
        Residue cys2 = new Residue("SER", 3);
        ala2.addAtom(new Atom("CA", Element.C, new double[] { 5.055, 74.864, 22.583 }));
        his4.addAtom(new Atom("CA", Element.C, new double[] { 7.321, 76.962, 20.326 }));
        cys2.addAtom(new Atom("CA", Element.C, new double[] { 6.020, 74.873, 17.386 }));
        container4 = GroupContainer.of(Stream.of(ala2, his4, cys2));
    }

    @Test
    public void shouldAlignFragments() throws IOException {
        // identical fragments with different orientation in the 3D space
        List<Protein> alignedFragments = Files.list(Paths.get(TestUtils.getResourceAsFilepath(FRAGMENT_DIR)))
                .map(ProteinParser::parsePDBFile)
                .limit(10)
                .collect(StructureCollectors.toAlignedEnsemble());

        // use first as reference
        Protein reference = alignedFragments.get(0);
        // all others should report a RMSD of 0.0 after alignment
        double maxRmsd = alignedFragments.stream()
                .skip(1)
                .mapToDouble(protein -> CoordinateUtils.calculateRMSD(reference, protein))
                .peek(System.out::println)
                .max()
                .orElse(Double.MAX_VALUE);

        Assert.assertEquals(0.0, maxRmsd, TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldAlignArbitraryPoints() {
        // compute alignment
        AlignmentResult alignmentResult = new SVDSuperimposer(AtomNameFilter.CA_ATOM_FILTER).align(container1, container2);
        System.out.println(Arrays.toString(alignmentResult.getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getRotation()));
        System.out.println("rmsd " + alignmentResult.getRmsd());
        Assert.assertEquals(0.19986, alignmentResult.getRmsd(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldAlignAnotherSetOfArbitraryPoints() {
        // compute alignment
        AlignmentResult alignmentResult = new SVDSuperimposer(AtomNameFilter.CA_ATOM_FILTER).align(container3, container4);
        System.out.println(Arrays.toString(alignmentResult.getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getRotation()));
        System.out.println("rmsd " + alignmentResult.getRmsd());
        Assert.assertEquals(0.0021, alignmentResult.getRmsd(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldSuperimposeBasedOnAlignment() {
        String initialCoordinates1 = container1.composePDBRecord();
        String initialCoordinates2 = container2.composePDBRecord();
        AlignmentResult alignmentResult = new SVDSuperimposer(AtomNameFilter.CA_ATOM_FILTER).align(container1, container2);
        double rmsd1 = alignmentResult.getRmsd();
        // coordinates should not be changed by aligning
        Assert.assertEquals(initialCoordinates1, container1.composePDBRecord());
        Assert.assertEquals(initialCoordinates2, container2.composePDBRecord());
        alignmentResult.transform(container2);
        double rmsd2 = CoordinateUtils.calculateRMSD(container1, container2);
        Assert.assertEquals(rmsd1, rmsd2, TestUtils.TOLERANT_ERROR_MARGIN);
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
        new SVDSuperimposer().align(protein1acj, protein1brr);
    }
}

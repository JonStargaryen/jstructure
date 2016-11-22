package alignment;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
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

import static util.TestUtils.FRAGMENT_DIR;

/**
 * Tests the capabilities of the SVDSuperimposer.
 * Created by S on 02.10.2016.
 */
public class SVDSuperimposerFunctionalTest {
    private Protein protein1acj;
    private Protein protein1brr;
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

        Group his1 = new Group("HIS", 1);
        Group asp1 = new Group("ASP", 2);
        Group ser1 = new Group("SER", 3);
        his1.addAtom(new Atom("CA", Element.C, new double[] { 6.994, 8.354, 42.405 }));
        asp1.addAtom(new Atom("CA", Element.C, new double[] { 9.429, 7.479, 48.266 }));
        ser1.addAtom(new Atom("CA", Element.C, new double[] { 5.547, 0.158, 42.050 }));
        container1 = Stream.of(his1, asp1, ser1).collect(StructureCollectors.toGroupContainer());

        Group his2 = new Group("HIS", 1);
        Group asp2 = new Group("ASP", 2);
        Group ser2 = new Group("SER", 3);
        his2.addAtom(new Atom("CA", Element.C, new double[] { 3.908, 12.066, -6.159 }));
        asp2.addAtom(new Atom("CA", Element.C, new double[] { 4.588, 6.531, -9.119 }));
        ser2.addAtom(new Atom("CA", Element.C, new double[] { 12.080, 12.645, -7.073 }));
        container2 = Stream.of(his2, asp2, ser2).collect(StructureCollectors.toGroupContainer());

        Group ala1 = new Group("HIS", 1);
        Group his3 = new Group("ASP", 2);
        Group cys1 = new Group("SER", 3);
        ala1.addAtom(new Atom("CA", Element.C, new double[] { 5.055, 74.865, 22.585 }));
        his3.addAtom(new Atom("CA", Element.C, new double[] { 7.320, 76.960, 20.325 }));
        cys1.addAtom(new Atom("CA", Element.C, new double[] { 6.021, 74.874, 17.385 }));
        container3 = Stream.of(ala1, his3, cys1).collect(StructureCollectors.toGroupContainer());

        Group ala2 = new Group("HIS", 1);
        Group his4 = new Group("ASP", 2);
        Group cys2 = new Group("SER", 3);
        ala2.addAtom(new Atom("CA", Element.C, new double[] { 5.055, 74.864, 22.583 }));
        his4.addAtom(new Atom("CA", Element.C, new double[] { 7.321, 76.962, 20.326 }));
        cys2.addAtom(new Atom("CA", Element.C, new double[] { 6.020, 74.873, 17.386 }));
        container4 = Stream.of(ala2, his4, cys2).collect(StructureCollectors.toGroupContainer());
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
        AlignmentResult alignmentResult = new SVDSuperimposer().align(container1, container2);
        System.out.println(Arrays.toString(alignmentResult.getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getRotation()));
        System.out.println("rmsd " + alignmentResult.getRmsd());
        Assert.assertEquals(0.19986, alignmentResult.getRmsd(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldAlignAnotherSetOfArbitraryPoints() {
        // compute alignment
        AlignmentResult alignmentResult = new SVDSuperimposer().align(container3, container4);
        System.out.println(Arrays.toString(alignmentResult.getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getRotation()));
        System.out.println("rmsd " + alignmentResult.getRmsd());
        Assert.assertEquals(0.0021, alignmentResult.getRmsd(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldSuperimposeBasedOnAlignment() {
        String initialCoordinates1 = container1.composePDBRecord();
        String initialCoordinates2 = container2.composePDBRecord();
        AlignmentResult alignmentResult = new SVDSuperimposer().align(container1, container2);
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

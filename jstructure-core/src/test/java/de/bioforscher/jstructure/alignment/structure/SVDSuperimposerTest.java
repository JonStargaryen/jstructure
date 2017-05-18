package de.bioforscher.jstructure.alignment.structure;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.parser.CIFParser;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Stream;

/**
 * Tests the capabilities of the SVDSuperimposer.
 * Created by S on 02.10.2016.
 */
public class SVDSuperimposerTest {
    private Protein protein1acj;
    private GroupContainer container1;
    private GroupContainer container2;
    private GroupContainer container3;
    private GroupContainer container4;

    private Group createGroup(String aminoAcidName, int residueNumber) {
        return new Group(aminoAcidName, residueNumber, CIFParser.parseLigandInformation(aminoAcidName), false);
    }

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParser.source("1acj").parse();

        Group his1 = createGroup("HIS", 1);
        Group asp1 = createGroup("ASP", 2);
        Group ser1 = createGroup("SER", 3);
        his1.addAtom(new Atom("CA", Element.C, new double[] { 6.994, 8.354, 42.405 }));
        asp1.addAtom(new Atom("CA", Element.C, new double[] { 9.429, 7.479, 48.266 }));
        ser1.addAtom(new Atom("CA", Element.C, new double[] { 5.547, 0.158, 42.050 }));
        container1 = Stream.of(his1, asp1, ser1).collect(StructureCollectors.toGroupContainer());

        Group his2 = createGroup("HIS", 1);
        Group asp2 = createGroup("ASP", 2);
        Group ser2 = createGroup("SER", 3);
        his2.addAtom(new Atom("CA", Element.C, new double[] { 3.908, 12.066, -6.159 }));
        asp2.addAtom(new Atom("CA", Element.C, new double[] { 4.588, 6.531, -9.119 }));
        ser2.addAtom(new Atom("CA", Element.C, new double[] { 12.080, 12.645, -7.073 }));
        container2 = Stream.of(his2, asp2, ser2).collect(StructureCollectors.toGroupContainer());

        Group ala3 = createGroup("HIS", 1);
        Group his3 = createGroup("ASP", 2);
        Group cys3 = createGroup("SER", 3);
        ala3.addAtom(new Atom("CA", Element.C, new double[] { 5.055, 74.865, 22.585 }));
        his3.addAtom(new Atom("CA", Element.C, new double[] { 7.320, 76.960, 20.325 }));
        cys3.addAtom(new Atom("CA", Element.C, new double[] { 6.021, 74.874, 17.385 }));
        container3 = Stream.of(ala3, his3, cys3).collect(StructureCollectors.toGroupContainer());

        Group ala4 = createGroup("HIS", 1);
        Group his4 = createGroup("ASP", 2);
        Group cys4 = createGroup("SER", 3);
        ala4.addAtom(new Atom("CA", Element.C, new double[] { 5.055, 74.864, 22.583 }));
        his4.addAtom(new Atom("CA", Element.C, new double[] { 7.321, 76.962, 20.326 }));
        cys4.addAtom(new Atom("CA", Element.C, new double[] { 6.020, 74.873, 17.386 }));
        container4 = Stream.of(ala4, his4, cys4).collect(StructureCollectors.toGroupContainer());
    }

    @Test
    public void shouldAlignArbitraryPoints() {
        // compute alignment
        StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(container1, container2);
        System.out.println(Arrays.toString(alignmentResult.getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getRotation()));
        System.out.println("rmsd " + alignmentResult.getAlignmentScore());
        System.out.println("rmsd2 " + LinearAlgebraAtom.calculateRmsd(alignmentResult.getOriginalReference(), alignmentResult.getOriginalQuery()));
        Assert.assertEquals(0.19986, alignmentResult.getAlignmentScore(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldAlignAnotherSetOfArbitraryPoints() {
        // compute alignment
        StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(container3, container4);
        System.out.println(Arrays.toString(alignmentResult.getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getRotation()));
        System.out.println("rmsd " + alignmentResult.getAlignmentScore());
        Assert.assertEquals(0.0021, alignmentResult.getAlignmentScore(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldSuperimposeBasedOnAlignment() {
        String initialCoordinates1 = container1.composePDBRecord();
        String initialCoordinates2 = container2.composePDBRecord();
        StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(container1, container2);
        double rmsd1 = alignmentResult.getAlignmentScore();
        // coordinates should not be changed by aligning
        Assert.assertEquals(initialCoordinates1, container1.composePDBRecord());
        Assert.assertEquals(initialCoordinates2, container2.composePDBRecord());
        alignmentResult.transform(container2);
        double rmsd2 = LinearAlgebraAtom.calculateRmsd(container1, container2);
        Assert.assertEquals(rmsd1, rmsd2, TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldResultInPerfectAlignment() {
        SVDSuperimposer svd = new SVDSuperimposer();
        StructureAlignmentResult result = svd.align(protein1acj, protein1acj);
        Assert.assertEquals(0.0, result.getAlignmentScore(), 0.001);
        Assert.assertArrayEquals(result.getTranslation(), AlignmentAlgorithm.NEUTRAL_TRANSLATION, 0.001);
    }

    @Test
    public void shouldResultInPerfectAlignmentForTransformedCopy() throws IOException {
        Protein protein1acjCopy = ProteinParser.source("1acj").parse();
        double[] translation = new double[] { 10, 20, 30 };
        LinearAlgebraAtom.transform(protein1acjCopy,
                translation,
                AlignmentAlgorithm.NEUTRAL_ROTATION);

        SVDSuperimposer svd = new SVDSuperimposer();
        StructureAlignmentResult result = svd.align(protein1acj.select()
                .aminoAcids()
                .asGroupContainer(), protein1acjCopy.select()
                .aminoAcids()
                .asGroupContainer());
        Assert.assertEquals(0.0, result.getAlignmentScore(), 0.001);
        Assert.assertArrayEquals(result.getTranslation(), LinearAlgebra3D.multiply(translation, -1.0), 0.001);
    }
}
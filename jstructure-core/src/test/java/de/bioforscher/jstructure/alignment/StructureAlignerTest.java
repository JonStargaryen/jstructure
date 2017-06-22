package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Stream;

/**
 * Tests for the structure aligner.
 * Created by bittrich on 6/19/17.
 */
public class StructureAlignerTest {
    private Protein protein1acj;
    private GroupContainer container1;
    private GroupContainer container2;
    private GroupContainer container3;
    private GroupContainer container4;

    private Group createGroup(String aminoAcidName, int residueNumber) {
        return new Group(aminoAcidName,
                new ResidueNumber(residueNumber),
                false);
    }

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParser.source("1acj").parse();

        Group his1 = createGroup("HIS", 1);
        Group asp1 = createGroup("ASP", 2);
        Group ser1 = createGroup("SER", 3);
        his1.addAtom(Atom.builder(Element.C, new double[] { 6.994, 8.354, 42.405 })
                .name("CA")
                .build());
        asp1.addAtom(Atom.builder(Element.C, new double[] { 9.429, 7.479, 48.266 })
                .name("CA")
                .build());
        ser1.addAtom(Atom.builder(Element.C, new double[] { 5.547, 0.158, 42.050 })
                .name("CA")
                .build());
        container1 = Stream.of(his1, asp1, ser1).collect(StructureCollectors.toGroupContainer());

        Group his2 = createGroup("HIS", 1);
        Group asp2 = createGroup("ASP", 2);
        Group ser2 = createGroup("SER", 3);
        his2.addAtom(Atom.builder(Element.C, new double[] { 3.908, 12.066, -6.159 })
                .name("CA")
                .build());
        asp2.addAtom(Atom.builder(Element.C, new double[] { 4.588, 6.531, -9.119 })
                .name("CA")
                .build());
        ser2.addAtom(Atom.builder(Element.C, new double[] { 12.080, 12.645, -7.073 })
                .name("CA")
                .build());
        container2 = Stream.of(his2, asp2, ser2).collect(StructureCollectors.toGroupContainer());

        Group ala3 = createGroup("HIS", 1);
        Group his3 = createGroup("ASP", 2);
        Group cys3 = createGroup("SER", 3);
        ala3.addAtom(Atom.builder(Element.C, new double[] { 5.055, 74.865, 22.585 })
                .name("CA")
                .build());
        his3.addAtom(Atom.builder(Element.C, new double[] { 7.320, 76.960, 20.325 })
                .name("CA")
                .build());
        cys3.addAtom(Atom.builder(Element.C, new double[] { 6.021, 74.874, 17.385 })
                .name("CA")
                .build());
        container3 = Stream.of(ala3, his3, cys3).collect(StructureCollectors.toGroupContainer());

        Group ala4 = createGroup("HIS", 1);
        Group his4 = createGroup("ASP", 2);
        Group cys4 = createGroup("SER", 3);
        ala4.addAtom(Atom.builder(Element.C, new double[] { 5.055, 74.864, 22.583 })
                .name("CA")
                .build());
        his4.addAtom(Atom.builder(Element.C, new double[] { 7.321, 76.962, 20.326 })
                .name("CA")
                .build());
        cys4.addAtom(Atom.builder(Element.C, new double[] { 6.020, 74.873, 17.386 })
                .name("CA")
                .build());
        container4 = Stream.of(ala4, his4, cys4).collect(StructureCollectors.toGroupContainer());
    }

    @Test
    public void shouldAlignArbitraryPoints() {
        // calculate alignment
        Alignment alignmentResult = StructureAligner.builder(container1, container2)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align();
        System.out.println(Arrays.toString(alignmentResult.getTransformation().getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getTransformation().getRotation()));
        System.out.println("rmsd " + alignmentResult.getAlignmentScore());
        Assert.assertEquals(0.19986, alignmentResult.getAlignmentScore(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldAlignAnotherSetOfArbitraryPoints() {
        // calculate alignment
        Alignment alignmentResult = StructureAligner.builder(container3, container4)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align();
        System.out.println(Arrays.toString(alignmentResult.getTransformation().getTranslation()));
        System.out.println(Arrays.deepToString(alignmentResult.getTransformation().getRotation()));
        System.out.println("rmsd " + alignmentResult.getAlignmentScore());
        Assert.assertEquals(0.0021, alignmentResult.getAlignmentScore(), TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void shouldManipulateCoordinatesInplace() {
        //TODO
    }

    @Test
    public void shouldNotManipulateCoordinatesCopy() {
        String initialCoordinates1 = container1.getPdbRepresentation();
        String initialCoordinates2 = container2.getPdbRepresentation();
        Alignment alignmentResult = StructureAligner.builder(container1, container2)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align();
        // coordinates should not be changed by aligning
        Assert.assertNotEquals(initialCoordinates2, alignmentResult.getAlignedQuery().getPdbRepresentation());
        Assert.assertEquals(initialCoordinates1, container1.getPdbRepresentation());
        Assert.assertEquals(initialCoordinates2, container2.getPdbRepresentation());
    }

    @Test
    public void shouldResultInPerfectAlignment() {
        Alignment alignmentResult = StructureAligner.builder(protein1acj, protein1acj)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align();
        Assert.assertEquals(0.0, alignmentResult.getAlignmentScore(), 0.001);
        Assert.assertArrayEquals(alignmentResult.getTransformation().getTranslation(), Transformation.NEUTRAL_TRANSLATION, 0.001);
    }

    @Test
    public void shouldResultInPerfectAlignmentForTransformedCopy() {
        Protein protein1acjCopy = ProteinParser.source("1acj").parse();
        double[] translation = new double[] { 10, 20, 30 };
        protein1acjCopy.calculate().transform(translation);

        Alignment alignmentResult = StructureAligner.builder(protein1acj, protein1acjCopy)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align();
        Assert.assertEquals(0.0, alignmentResult.getAlignmentScore(), 0.001);
        Assert.assertArrayEquals(alignmentResult.getTransformation().getTranslation(), LinearAlgebra.on(translation).multiply(-1.0).getValue(), 0.001);
    }

    @Test
    public void shouldAlignSyntheticContainers() {
        double rmsd = StructureAligner.builder(protein1acj.select()
                .aminoAcids()
                .asGroupContainer(), protein1acj.select()
                .aminoAcids()
                .asGroupContainer())
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align()
                .getAlignmentScore();
        Assert.assertEquals(0.0, rmsd, TestUtils.TOLERANT_ERROR_MARGIN);
    }
}
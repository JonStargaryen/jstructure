package de.bioforscher.jstructure.parser;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.junit.Assert;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static org.hamcrest.CoreMatchers.endsWith;
import static org.hamcrest.CoreMatchers.startsWith;

/**
 * Checks functions of the PDB parser and its integrity with the data model.
 * Created by S on 29.09.2016.
 */
public class ProteinParserTest {
    private static final String PDB_EXTENSION = ".pdb";
    private static final String PDB_DIRECTORY = "parser/";
    private static final List<String> PDB_IDS = Arrays.asList("2w0l", "1acj", "1asz", "4bpm");
    private static final List<Path> PDB_FILE_PATHS = PDB_IDS.stream()
            .map(f -> PDB_DIRECTORY + f + PDB_EXTENSION)
            .map(TestUtils::getResourceAsFilepath)
            .map(Paths::get)
            .collect(Collectors.toList());
    /**
     * contains selenomethionine at pos 1, marked as HETATM
     */
    private static final String NON_STANDARD_PDB_ID = "1dw9";

    @Test
    public void shouldHandleProteinWithOligopeptide() {
        Protein protein = ProteinParser.source("1lyb").parse();
        // contains peptide-like groups - ought to be an AminoAcid
        AminoAcid statine = (AminoAcid) protein.select()
                .groupName("STA")
                .asGroup();
        Assert.assertTrue(statine.getPolymerType() == GroupPrototype.PolymerType.PEPTIDE_LIKE);
        Assert.assertTrue(statine.isAminoAcid());
        Assert.assertFalse(statine.isLigand());
        Assert.assertTrue(protein.aminoAcids().collect(Collectors.toList()).contains(statine));
    }

    @Test
    public void shouldParseNonStandardAminoAcid() {
        ProteinParser.source(Paths.get(TestUtils.getResourceAsFilepath(PDB_DIRECTORY + "nonstandard/1dw9-first-selenomethionine.pdb"))).parse();
    }

    @Test
    public void shouldParseInsertedAminoAcids() {
        Protein protein = ProteinParser.source("2w0l").parse();
        List<Group> groups = protein.select()
                .chainName("A")
                .residueNumber(95)
                .asFilteredGroups()
                .collect(Collectors.toList());
        Assert.assertEquals(2, groups.size());
        groups.forEach(System.out::println);
    }

    @Test
    public void shouldParseAlternativePositions() {
        Protein protein = ProteinParser.source("4bpm").parse();
        protein.select()
                // atom 61 an alternative position
                .pdbSerial(61)
                .asAtom();
    }

    @Test
    public void shouldHandleProteinWithNonStandardAminoAcids() {
        Protein protein = ProteinParser.source(NON_STANDARD_PDB_ID).parse();
        // ensure that the initial selenomethionine stored as HETATM is correctly parsed
        System.out.println(protein.getAminoAcidSequence());
        Assert.assertThat(protein.getAminoAcidSequence(), startsWith("M"));
    }

    /**
     * Tests whether the parser can handle some usual structures.
     */
    @Test
    public void shouldParseAllStructureFiles() {
        PDB_FILE_PATHS.stream()
                .map(ProteinParser::source)
                .forEach(ProteinParser.OptionalSteps::parse);
    }

    @Test
    public void shouldFetchPdbStructureFromPDB() {
        PDB_IDS.forEach(id -> ProteinParser.source(id).parse());
    }

    @Test(expected = ParsingException.class)
    public void shouldFailForInvalidStructure() {
        ProteinParser.source(Paths.get(TestUtils.getResourceAsFilepath("pdb/invalid.pdb"))).strictMode(true).parse();
    }

    @Test
    public void shouldHandleModifiedResidue() {
        Protein protein = ProteinParser.source("1brr").parse();

        Group pca = protein.select()
                .chainName("C")
                .residueNumber(1)
                .asGroup();
        Assert.assertTrue(pca.isAminoAcid());
        Assert.assertFalse(pca.isLigand());
        // assert correct mapping of PCA to GLU
        Assert.assertEquals("incorrect mapping of PCA to GLU",
                "E",
                pca.getGroupPrototype().getOneLetterCode().get());
    }

    @Test
    public void shouldAnnotateHetAtmsCorrectlyFor1bs2() {
        /*
         * 1bs2 is an aars structure with the amino acid arginine in the binding site (annotated as ATOM record), some
         * water (annotated as HETATM)
         */
        Protein protein1bs2 = ProteinParser.source("1bs2").parse();

        List<Group> waters = protein1bs2.select()
                .water()
                .asFilteredGroups()
                .collect(Collectors.toList());

        waters.forEach(group -> {
            Assert.assertTrue(group.isLigand());
            Assert.assertTrue("water records ought to start with HETATM", group.getPdbRepresentation().startsWith(Atom.HETATM_PREFIX));
        });

        Group arginineAsLigand = protein1bs2.select()
                // in chain A at resNum 900 there is the ARG ligand
                .residueNumber(900)
                .asGroup();

        // assert that selection does not return ARG ligand as normal amino acid
        boolean arginineLigandIsNoAminoAcid = protein1bs2.aminoAcids()
                .noneMatch(group -> group.equals(arginineAsLigand));
        Assert.assertTrue("amino acid ligand ought to be not a part of the amino acid chain", arginineLigandIsNoAminoAcid);

        // ensure last amino acid is MET and not the ARG ligand
        Assert.assertThat(protein1bs2.getAminoAcidSequence(), endsWith("M"));

        List<Group> hetatm1bs2 = protein1bs2.select()
                .hetatms()
                .asFilteredGroups()
                .collect(Collectors.toList());

        Assert.assertTrue(hetatm1bs2.containsAll(waters) && hetatm1bs2.contains(arginineAsLigand));
    }

    /**
     * Tests whether <tt>ATOM</tt> records are parsed and written correctly.
     */
    @Test
    public void shouldWriteEqualAtomRecords() {
        PDB_IDS.forEach(this::checkAgreement);
    }

    private void checkAgreement(String pdbId) {
        System.out.println("checking agreement between written and expected ATOM records for " + pdbId);
        Protein protein = ProteinParser.source(pdbId).parse();
        List<String> writtenLines = Pattern.compile("\n")
                .splitAsStream(protein.getPdbRepresentation())
                .collect(Collectors.toList());
        try {
            List<String> expectedLines = Files.lines(Paths.get(TestUtils.getResourceAsFilepath("parser/parsed/" + pdbId + ".pdb")))
                    .collect(Collectors.toList());

            for(int i = 0; i < writtenLines.size(); i++) {
                // some file have shorter lines
                if(expectedLines.get(i).length() > writtenLines.get(i).length()) {
                    Assert.assertTrue("ATOM records do not match!" + System.lineSeparator() +
                            "expected: " + expectedLines.get(i) + System.lineSeparator() +
                            "actual:   " + writtenLines.get(i), expectedLines.get(i).startsWith(writtenLines.get(i)));
                } else {
                    Assert.assertTrue("ATOM records do not match!" + System.lineSeparator() +
                            "expected: " + expectedLines.get(i) + System.lineSeparator() +
                            "actual:   " + writtenLines.get(i), writtenLines.get(i).startsWith(expectedLines.get(i)));
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException("did not find expected output file for " + pdbId, e);
        }
    }
}

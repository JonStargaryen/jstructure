package parser;

import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ParsingException;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Test;
import util.TestUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.hamcrest.CoreMatchers.endsWith;
import static org.hamcrest.CoreMatchers.startsWith;

/**
 * Checks functions of the PDB parser and its integrity with the data de.bioforscher.explorer.helices.model.
 * Created by S on 29.09.2016.
 */
public class ProteinParserFunctionalTest {
    public static final String PDB_EXTENSION = ".pdb";
    public static final String PDB_DIRECTORY = "parser/";
    private static final String PDB_PARSED_DIRECTORY = PDB_DIRECTORY + "parsed/";
    private static final List<String> PDB_IDS = Arrays.asList("1acj", "1asz", "1brr", "4cha");
    private static final List<String> PDB_FILE_PATHS = PDB_IDS.stream()
            .map(f -> PDB_DIRECTORY + f + PDB_EXTENSION).collect(Collectors.toList());
    /**
     * contains selenomethionine at pos 1, marked as HETATM
     */
    private static final String NON_STANDARD_PDB_ID = "1dw9";

    @Test
    public void shouldParseNonStandardAminoAcid() {
        parseFilepath(PDB_DIRECTORY + "nonstandard/1dw9-first-selenomethionine.pdb");
    }

    @Test
    public void shouldHandleProteinWithNonStandardAminoAcids() {
        Protein protein = ProteinParser.parseProteinById(NON_STANDARD_PDB_ID);
        // ensure that the initial selenomethionine stored as HETATM is correctly parsed
        System.out.println(protein.getAminoAcidSequence());
        Assert.assertThat(protein.getAminoAcidSequence(), startsWith("M"));
    }

    /**
     * Tests whether the parser can handle some usual structures.
     */
    @Test
    public void shouldParseAllStructureFiles() {
        PDB_FILE_PATHS.forEach(ProteinParserFunctionalTest::parseFilepath);
    }

    @Test
    public void shouldFetchPdbStructureFromPDB() {
        PDB_IDS.forEach(ProteinParser::parseProteinById);
    }

    @Test(expected = ParsingException.class)
    public void shouldFailForInvalidStructure() {
        parseFilepath("pdb/invalid.pdb");
    }

    @Test
    public void shouldAnnotateHetAtmsCorrectlyFor1bs2() {
        /*
         * 1bs2 is an aars structure with the amino acid arginine in the binding site (annotated as ATOM record), some
         * water (annotated as HETATM)
         */
        Protein protein1bs2 = ProteinParser.parseProteinById("1bs2");

        List<Group> waters = Selection.on(protein1bs2)
                .water()
                .asFilteredGroups()
                .collect(Collectors.toList());

        waters.forEach(group -> {
            Assert.assertTrue(group.isLigand());
            Assert.assertTrue(group.composePDBRecord().startsWith(ProteinParser.HETATM_PREFIX));
        });

        Group arginineAsLigand = Selection.on(protein1bs2)
                // in chain A at resNum 900 there is the ARG ligand
                .residueNumber(900)
                .asGroup();

        // assert that selection does not return ARG ligand as normal amino acid
        boolean arginineLigandIsNoAminoAcid = Selection.on(protein1bs2)
                .aminoAcids()
                .asFilteredGroups()
                .noneMatch(group -> group.equals(arginineAsLigand));
        Assert.assertTrue(arginineLigandIsNoAminoAcid);

        // ensure last amino acid is MET and not the ARG ligand
        Assert.assertThat(protein1bs2.getAminoAcidSequence(), endsWith("M"));

        List<Group> hetatm1bs2 = Selection.on(protein1bs2)
                .hetatms()
                .asFilteredGroups()
                .collect(Collectors.toList());

        Assert.assertTrue(hetatm1bs2.containsAll(waters) && hetatm1bs2.contains(arginineAsLigand));
    }

    /**
     * Tests whether <tt>ATOM</tt> records findAny parsed and written correctly.
     */
    @Test
    public void shouldWriteEqualAtomRecords() {
        for (String pdbFilePath : PDB_FILE_PATHS) {
            Protein protein = parseFilepath(pdbFilePath);
            protein.chains().forEach(chain -> System.out.println("parsed getChain '" + chain.getChainId() + "'"));
            String parsedFilepath = TestUtils.getResourceAsFilepath(pdbFilePath.replace(PDB_DIRECTORY, PDB_PARSED_DIRECTORY));
            try {
                List<String> originalLines = Files.readAllLines(new File(parsedFilepath).toPath());
                List<String> writtenLines = Arrays.asList(protein.composePDBRecord().split(System.lineSeparator()));
                if(originalLines.size() != writtenLines.size()) {
                    System.out.println("original and written PDB file differ in size!");
                    originalLines.removeAll(writtenLines);
                    System.out.println("missing from original: " + originalLines);
                    Assert.fail();
                }
                for(int i = 0; i < originalLines.size(); i++) {
                    String originalLine = originalLines.get(i);
                    String writtenLine = writtenLines.get(i);
                    if(originalLine.length() != writtenLine.length()) {
                        System.out.println("lines differ in size!");
                        System.out.println("'" + originalLine + "'");
                        System.out.println("'" + writtenLine + "'");
                        Assert.fail();
                    }
                    Assert.assertEquals(originalLine, writtenLine);
                }
            } catch (IOException e) {
                Assert.fail("failed with IOException: " + e.getLocalizedMessage());
            }
        }
    }

    public static Protein parseFilepath(String filepath) {
        System.out.println("parsing PDB file " + filepath);
        return ProteinParser.parsePDBFile(TestUtils.getResourceAsFilepath(filepath));
    }
}

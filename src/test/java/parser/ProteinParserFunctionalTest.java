package parser;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import util.TestUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by S on 29.09.2016.
 */
public class ProteinParserFunctionalTest {
    private static final String PDB_EXTENSION = ".pdb";
    private static final String PDB_DIRECTORY = "pdb/";
    private static final String PDB_PARSED_DIRECTORY = PDB_DIRECTORY + "parsed/";
    private static final List<String> PDB_IDS = Arrays.asList("1acj", "1asz", "1brr", "4cha");
    private static final List<String> PDB_FILE_PATHS = PDB_IDS.stream()
            .map(f -> PDB_DIRECTORY + f + PDB_EXTENSION).collect(Collectors.toList());
    private ProteinParser proteinParser;

    @Before
    public void setup() {
        this.proteinParser = new ProteinParser();
    }

    /**
     * Tests whether the parser can handle some usual structures.
     */
    @Test
    public void shouldParseAllStructureFiles() {
        PDB_FILE_PATHS.forEach(this::parseFilepath);
    }

    @Test
    public void shouldFetchPdbStructureFromPDB() {
        PDB_IDS.forEach(i -> {
            try {
                proteinParser.parseProteinById(i);
            } catch (IOException e) {
                Assert.fail("failed with IOException: " + e.getLocalizedMessage());
            }
        });
    }

    /**
     * Tests whether <tt>ATOM</tt> records get parsed and written correctly.
     */
    @Test
    public void shouldWriteEqualAtomRecords() {
        for (String pdbFilePath : PDB_FILE_PATHS) {
            Protein protein = parseFilepath(pdbFilePath);
            protein.chains().forEach(chain -> System.out.println("parsed chain '" + chain.getChainId() + "'"));
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

    private Protein parseFilepath(String filepath) {
        System.out.println("parsing PDB file " + filepath);
        try {
            return this.proteinParser.parsePDBFile(TestUtils.getResourceAsFilepath(filepath));
        } catch (IOException e) {
            Assert.fail("failed with IOException: " + e.getLocalizedMessage());
        }
        // compiler please...
        return null;
    }
}
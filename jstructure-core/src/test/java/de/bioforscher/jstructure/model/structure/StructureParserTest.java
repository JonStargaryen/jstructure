package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.nucleotide.Nucleotide;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.hamcrest.CoreMatchers.endsWith;
import static org.hamcrest.CoreMatchers.startsWith;

/**
 * Checks functions of the PDB parser and its integrity with the data model.
 * Created by S on 29.09.2016.
 */
public class StructureParserTest {
    private static final String PDB_DIRECTORY = "parser/";
    private static final List<TestUtils.SupportedProtein> PDB_IDS = Stream.of(TestUtils.SupportedProtein.values())
        .collect(Collectors.toList());
    /**
     * contains selenomethionine at pos 1, marked as HETATM
     */
    private static final String NON_STANDARD_PDB_ID = "1dw9";

    @Test
    public void shouldHandleRNA() {
        String pdbId = "5new";
        Structure structure = StructureParser.fromPdbId(pdbId)
                .minimalParsing(true)
                .parse();
        Assert.assertEquals("number of nucleotide groups did not match", 8, structure.nucleotides().count());
        Assert.assertEquals("number of atoms did not match", 1261, structure.atoms().count());
        structure.nucleotides()
                .map(Nucleotide::getN1)
                .forEach(Assert::assertNotNull);
    }

    @Test
    public void shouldHandleLegacyHeaders() {
        String pdbId = "2m7g";
        Structure structure = StructureParser.fromPdbId(pdbId)
                .minimalParsing(true)
                .parse();
        Assert.assertEquals("pdbId does not match expectation", pdbId, structure.getProteinIdentifier().getPdbId());
        Assert.assertEquals("date does not match", "2013-04-22", structure.getDepositionDate().toString());
    }

    @Test
    public void shouldFormatPdbHeaderCorrectly() {
        Structure protein = StructureParser.fromInputStream(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1AR1))
                .minimalParsing(true)
                .parse();
        String expected = "HEADER    COMPLEX (OXIDOREDUCTASE/ANTIBODY)       08-AUG-97   1AR1              " + System.lineSeparator();
        String actual = protein.getHeader();
        Assert.assertEquals("length of header do not match", expected.length(), actual.length());
        Assert.assertEquals("header line did not match expectation", expected, actual);
    }

    @Test
    public void shouldSkipLigandParsing() {
        String id = "5GRO";
        Structure proteinFast = StructureParser.fromPdbId(id).minimalParsing(true).parse();
        List<Group> ligandsFast = proteinFast.select()
                .ligands()
                .asFilteredGroups()
                .collect(Collectors.toList());

        Structure proteinConventional = StructureParser.fromPdbId(id).parse();
        List<Group> ligandsConventional = proteinConventional.select()
                .ligands()
                .asFilteredGroups()
                .collect(Collectors.toList());

        Assert.assertEquals("amino acids should match in size", proteinConventional.aminoAcids().count(), proteinFast.aminoAcids().count());
        Assert.assertEquals("nucleotides should match in size", proteinConventional.select().nucleotides().asFilteredGroups().count(), proteinFast.select().nucleotides().asFilteredGroups().count());
        Assert.assertEquals("annotated ligands should match", ligandsConventional.size(), ligandsFast.size());
    }

    @Test
    public void shouldHandleProteinWithOligopeptide() {
        Structure protein = StructureParser.fromPdbId("1lyb").parse();
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
        StructureParser.fromInputStream(TestUtils.getResourceAsInputStream(PDB_DIRECTORY + "nonstandard/1dw9-first-selenomethionine.pdb")).parse();
    }

    @Test
    public void shouldRetrieveAminoAcidSequenceOfAllChains() {
        Structure protein = StructureParser.fromInputStream(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_5OAZ))
                .minimalParsing(true)
                .parse();
        String sequence = protein.getAminoAcidSequence();
        long numberOfGroups = protein.chains()
                .flatMap(Chain::groups)
                .filter(Group::isAminoAcid)
                .count();
        Assert.assertEquals("sequence string is not of expected size", numberOfGroups, sequence.length());
    }

    @Test
    public void shouldParseInsertedAminoAcids() {
        Structure protein = StructureParser.fromPdbId("2w0l").parse();
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
        Structure protein = StructureParser.fromPdbId("4bpm").parse();
        protein.select()
                // atom 61 an alternative position
                .pdbSerial(61)
                .asAtom();
    }

    @Test
    public void shouldHandleProteinWithNonStandardAminoAcids() {
        Structure protein = StructureParser.fromPdbId(NON_STANDARD_PDB_ID).parse();
        // ensure that the initial selenomethionine stored as HETATM is correctly parsed
        System.out.println(protein.getAminoAcidSequence());
        Assert.assertThat(protein.getAminoAcidSequence(), startsWith("M"));
    }

    @Test(expected = ParsingException.class)
    public void shouldFailForInvalidStructure() {
        StructureParser.fromInputStream(TestUtils.getResourceAsInputStream("pdb/invalid.pdb")).strictMode(true).parse();
    }

    @Test
    public void shouldHandleModifiedResidue() {
        Structure protein = StructureParser.fromPdbId("1brr").parse();

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
        Structure protein1bs2 = StructureParser.fromPdbId("1bs2").parse();

        List<Group> waters = protein1bs2.select()
                .water()
                .asFilteredGroups()
                .collect(Collectors.toList());

        waters.forEach(group -> {
            Assert.assertTrue(group.isLigand());
            Assert.assertTrue("water records ought to start with HETATM",
                    group.getPdbRepresentation().contains(Atom.HETATM_PREFIX));
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
        //TODO remove and add other test-cases
//        PDB_IDS.forEach(this::checkAgreement);
        checkAgreement(TestUtils.SupportedProtein.PDB_3G1H);
    }

    private void checkAgreement(TestUtils.SupportedProtein supportedProtein) {
        String pdbId = supportedProtein.name().split("_")[1];

        System.out.println("checking agreement between written and expected ATOM records for " + pdbId);
        Structure protein = StructureParser.fromInputStream(TestUtils.getProteinInputStream(supportedProtein)).parse();
        List<String> writtenLines = Pattern.compile("\n")
                .splitAsStream(protein.getPdbRepresentation())
                .collect(Collectors.toList());
        List<String> expectedLines = TestUtils.getResourceAsLines("parser/parsed/" + pdbId.toLowerCase() + ".pdb");

        Assert.assertEquals("number of lines do not match", expectedLines.size(), writtenLines.size());

        for(int i = 0; i < writtenLines.size(); i++) {
            System.out.println(writtenLines.get(i));
            // some files have shorter lines
            if(expectedLines.get(i).length() > writtenLines.get(i).length()) {
                Assert.assertTrue("ATOM records do not match" + System.lineSeparator() +
                        "expected: " + expectedLines.get(i) + System.lineSeparator() +
                        "actual:   " + writtenLines.get(i), expectedLines.get(i).startsWith(writtenLines.get(i)));
            } else {
                Assert.assertTrue("ATOM records do not match" + System.lineSeparator() +
                        "expected: " + expectedLines.get(i) + System.lineSeparator() +
                        "actual:   " + writtenLines.get(i), writtenLines.get(i).startsWith(expectedLines.get(i)));
            }
        }
    }

    @Test
    public void shouldHandleLargeCoordinatesFor3hqv() {
        Structure structure = StructureParser.fromPdbId("3hqv").parse();
        Pattern.compile("\n").splitAsStream(structure.getPdbRepresentation())
                .filter(line -> line.startsWith("ATOM") || line.startsWith("HETATM"))
                .forEach(line -> Assert.assertFalse("file formatted wrong - contains ',' - line: " + line,
                        line.contains(",")));
    }

    @Test
    public void shouldCorrectOffsetInHeaderFor4xdb() {
        Structure structure = StructureParser.fromPdbId("4xdb").parse();
        String header = structure.getHeader();

        Assert.assertEquals("header not fixed",
                "HEADER    OXIDOREDUCTASE, MEMBRANE PROTEIN, FLAVOP19-DEC-14   4XDB              ".trim(),
                header.trim());
    }

    @Test
    public void shouldLeaveCorrectHeaderUntouchedFor1acj() {
        Structure structure = StructureParser.fromPdbId("1acj").parse();
        String header = structure.getHeader();

        Assert.assertEquals("header was manipulated",
                "HEADER    HYDROLASE(CARBOXYLIC ESTERASE)          18-AUG-93   1ACJ              ".trim(),
                header.trim());
    }
}

package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.LigandContactScreener;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test for the ligand screener.
 * Created by bittrich on 7/13/17.
 */
public class LigandContactScreenerImplTest {
    private LigandContactScreener ligandContactScreener;
    private Structure protein;
    private AminoAcid aminoAcidNextToLigand;
    private AminoAcid aminoAcidNextToWater;

    @Before
    public void setup() {
        ligandContactScreener = new LigandContactScreenerImpl();
        protein = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        aminoAcidNextToLigand = protein.select()
                .chainName("A")
                .residueNumber(84)
                .asAminoAcid();
        aminoAcidNextToWater = protein.select()
                .chainName("A")
                .residueNumber(416)
                .asAminoAcid();
    }

    @Test
    public void shouldFindNeighboredLigands() {
        Assert.assertTrue("did not find any ligand contacts",
                ligandContactScreener.determineNumberOfLigandContacts(protein, aminoAcidNextToLigand) > 0);
    }

    @Test
    public void shouldIgnoreWater() {
        // there is water at 416 but no actual ligand
        Assert.assertTrue("probably considered water as ligand",
                ligandContactScreener.determineNumberOfLigandContacts(protein, aminoAcidNextToWater) == 0);
    }
}
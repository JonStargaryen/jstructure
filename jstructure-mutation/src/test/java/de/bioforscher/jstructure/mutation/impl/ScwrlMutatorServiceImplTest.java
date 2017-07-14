package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Test for the scwrl integration.
 * Created by bittrich on 7/14/17.
 */
public class ScwrlMutatorServiceImplTest {
    private ScwrlMutatorServiceImpl scwrlMutatorService;
    private String pdbId;
    private Protein protein;
    private String chainId;
    private int residueNumber;
    private AminoAcid aminoAcidToMutate;
    private AminoAcid.Family targetAminoAcid;

    @Before
    public void setup() {
        scwrlMutatorService = new ScwrlMutatorServiceImpl();
        pdbId = "5oaz";
        protein = ProteinParser.source(pdbId)
                .minimalParsing(true)
                .parse();
        chainId = "A";
        residueNumber = 80;
        aminoAcidToMutate = protein.select()
                .chainName(chainId)
                .residueNumber(residueNumber)
                .asAminoAcid();
        targetAminoAcid = AminoAcid.Family.ALANINE;
    }

    @Test
    @Ignore
    public void shouldIntroduceMutation() throws IOException {
        Protein mutatedProtein = scwrlMutatorService.mutateAminoAcid(protein, aminoAcidToMutate, targetAminoAcid);
        //TODO global test dirs
        Files.write(Paths.get("/home/bittrich/Downloads/original.pdb"), protein.getPdbRepresentation().getBytes());
        Files.write(Paths.get("/home/bittrich/Downloads/mutated.pdb"), mutatedProtein.getPdbRepresentation().getBytes());
    }

    @Test
    public void shouldProvideSequenceString() {
        String originalSequence = protein.getAminoAcidSequence();
        String mutatedSequence = scwrlMutatorService.composeMutateScwrlSequence(protein, aminoAcidToMutate, targetAminoAcid);
        System.out.println(originalSequence);
        System.out.println(mutatedSequence);
        Assert.assertEquals("swrcl sequence did not match expectation", "dpnlfvalydfvasgdntAsitkgeklrvlgynhngewceaqtkngqgwvpsnyitpvnlfvalydfvasgdntlsitkgeklrvlgynhngewceaqtkngqgwvpsnyitpvn", mutatedSequence);
    }
}
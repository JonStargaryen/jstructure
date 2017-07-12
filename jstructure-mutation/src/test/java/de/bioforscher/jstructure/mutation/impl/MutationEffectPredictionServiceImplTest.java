package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.mutation.MutationJob;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.PrimaryUniProtAccession;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

/**
 * Test for the function of the mutation effect predictor.
 * Created by bittrich on 7/11/17.
 */
public class MutationEffectPredictionServiceImplTest {
    private MutationEffectPredictionServiceImpl mutationEffectPredictor;

    @Before
    public void setup() {
        mutationEffectPredictor = new MutationEffectPredictionServiceImpl();
    }

    @Test
    public void shouldBlastForHomologousSequenceAndStructures() throws ExecutionException {
        MutationJob mutationJob = new MutationJobImpl("CYC32_DESDN",
                "ETFEIPESVTMSPKQFEGYTPKKGDVTFNHASHMDIACQQCHHTVPDTYTIESCMTEGCHDNIKERTEISSVYRTFHTTKDSEKSCVGCHRELKRQGPSDAPLACNSCHVQ");
        mutationEffectPredictor.createMultiSequenceAlignment(mutationJob);

        System.out.println("homologous sequences:");
        System.out.println(mutationJob.getHomologousEntryContainer()
                .getUniProtEntries()
                .stream()
                .map(UniProtEntry::getPrimaryUniProtAccession)
                .map(PrimaryUniProtAccession::getValue)
                .collect(Collectors.toList()));
        Assert.assertTrue(mutationJob.getHomologousEntryContainer().getUniProtEntries().size() >= 8);

        System.out.println("homologous pdb chains:");
        System.out.println(mutationJob.getHomologousPdbChains()
                .stream()
                .map(Chain::getChainIdentifier)
                .map(ChainIdentifier::getFullName)
                .collect(Collectors.toList()));
        Assert.assertTrue(mutationJob.getHomologousPdbChains().size() > 0);
    }

    @Test
    public void shouldExtractFragmentsAndRunItemsetMiner() throws ExecutionException {
        MutationJob mutationJob = new MutationJobImpl("CYC32_DESDN",
                "ETFEIPESVTMSPKQFEGYTPKKGDVTFNHASHMDIACQQCHHTVPDTYTIESCMTEGCHDNIKERTEISSVYRTFHTTKDSEKSCVGCHRELKRQGPSDAPLACNSCHVQ");
        mutationEffectPredictor.createMultiSequenceAlignment(mutationJob);
        System.out.println(mutationJob.composeMutationDescriptor(73, AminoAcid.Family.ALANINE));
    }
}
package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.mutation.MutationEffectPrediction;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.PrimaryUniProtAccession;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

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
    public void shouldBlastForHomologousSequenceAndStructures() {
        MutationEffectPrediction mutationEffectPrediction = new MutationEffectPredictionImpl("CYC32_DESDN",
                "ETFEIPESVTMSPKQFEGYTPKKGDVTFNHASHMDIACQQCHHTVPDTYTIESCMTEGCHDNIKERTEISSVYRTFHTTKDSEKSCVGCHRELKRQGPSDAPLACNSCHVQ");
        mutationEffectPredictor.executeBlastQuery(mutationEffectPrediction);

        System.out.println("homologous sequences:");
        System.out.println(mutationEffectPrediction.getHomologousEntryContainer()
                .getUniProtEntries()
                .stream()
                .map(UniProtEntry::getPrimaryUniProtAccession)
                .map(PrimaryUniProtAccession::getValue)
                .collect(Collectors.toList()));
        Assert.assertTrue(mutationEffectPrediction.getHomologousEntryContainer().getUniProtEntries().size() >= 8);

        System.out.println("homologous pdb chains:");
        System.out.println(mutationEffectPrediction.getHomologousPdbChains()
                .stream()
                .map(Chain::getChainIdentifier)
                .map(ChainIdentifier::getFullName)
                .collect(Collectors.toList()));
        Assert.assertTrue(mutationEffectPrediction.getHomologousPdbChains().size() > 0);
    }

    @Test
    public void shouldExtractFragmentsAndRunItemsetMiner() {
        MutationEffectPrediction mutationEffectPrediction = new MutationEffectPredictionImpl("CYC32_DESDN",
                "ETFEIPESVTMSPKQFEGYTPKKGDVTFNHASHMDIACQQCHHTVPDTYTIESCMTEGCHDNIKERTEISSVYRTFHTTKDSEKSCVGCHRELKRQGPSDAPLACNSCHVQ");
        mutationEffectPredictor.executeBlastQuery(mutationEffectPrediction);
        mutationEffectPrediction.predictMutationEffect(50, AminoAcid.Family.ALANINE);
    }
}
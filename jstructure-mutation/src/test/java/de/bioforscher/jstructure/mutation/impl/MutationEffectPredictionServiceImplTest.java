package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.mutation.MutationEffectPredictionService;
import de.bioforscher.jstructure.mutation.MutationJob;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test for the mutation effect predictor.
 * Created by bittrich on 7/14/17.
 */
public class MutationEffectPredictionServiceImplTest {
    private MutationEffectPredictionService mutationEffectPredictionService;
    private Chain referenceChain;

    @Before
    public void setup() {
        mutationEffectPredictionService = new MutationEffectPredictionServiceImpl();
        referenceChain = ProteinParser.source("1BTA")
                .minimalParsing(true)
                .parse()
                .select()
                .chainName("A")
                .asChain();
    }

    @Test
    public void shouldCreateMutationJob() {
        MutationJob mutationJob = mutationEffectPredictionService.createMutationJob("test-job", referenceChain);
        Assert.assertNotNull(mutationJob);
    }
}
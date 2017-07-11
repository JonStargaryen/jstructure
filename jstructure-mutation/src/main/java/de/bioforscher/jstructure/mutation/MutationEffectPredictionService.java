package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.mutation.impl.MutationEffectPredictionServiceImpl;

/**
 * The wrapping class which implements the workflow to predict the effect of mutations.
 * Created by bittrich on 7/11/17.
 */
public interface MutationEffectPredictionService {
    MutationEffectPredictionService INSTANCE = new MutationEffectPredictionServiceImpl();

    /**
     * Create the object to predict the effect of a mutation for a given sequence.
     * @param identifier this job's name
     * @param sequence the sequence to process
     */
    MutationEffectPrediction getMutationEffectPrediction(String identifier, String sequence);
}

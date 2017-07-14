package de.bioforscher.jstructure.mutation.old;

import java.util.concurrent.ExecutionException;

/**
 * The wrapping class which implements the workflow to predict the effect of mutations.
 * Created by bittrich on 7/11/17.
 */
@Deprecated
public interface MutationEffectPredictionService {
    /**
     * Create the object to predict the effect of a mutation for a given sequence.
     * @param identifier this job's name
     * @param sequence the sequence to process
     * @throws ExecutionException when an internal task struggles to a degree compromising this task's integrity
     */
    MutationJob createMutationJob(String identifier, String sequence) throws ExecutionException;
}

package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.mutation.impl.SimpleMutatorServiceImpl;

/**
 * Creates the data structure to predict the effect of a mutation.
 * Created by bittrich on 7/13/17.
 */
public class MutationFeatureVectorFactory {
    private static final MutationFeatureVectorFactory INSTANCE = new MutationFeatureVectorFactory();
    private final MutatorService mutatorService;

    private MutationFeatureVectorFactory() {
        this.mutatorService = new SimpleMutatorServiceImpl();
    }

    public static MutationFeatureVectorFactory getInstance() {
        return INSTANCE;
    }
}

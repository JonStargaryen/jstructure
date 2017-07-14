package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Predicts the effect of a given mutation in the context of a protein structure.
 * Created by bittrich on 7/13/17.
 */
public interface MutationEffectPredictionService {
    MutationJob createMutationJob(String jobName, Chain referenceChain);

    MutationFeatureVector createMutationFeatureVector(MutationJob mutationJob,
                                                      ResidueIdentifier residueIdentifierToMutate,
                                                      AminoAcid.Family mutationTarget);
}

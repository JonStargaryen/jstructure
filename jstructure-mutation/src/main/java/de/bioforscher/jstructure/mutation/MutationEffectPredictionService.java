package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

/**
 * Predicts the effect of a given mutation in the context of a protein structure.
 * Created by bittrich on 7/13/17.
 */
public interface MutationEffectPredictionService {
    MutationJob createMutationJob(String jobName,
                                  Chain referenceChain,
                                  LocalBlastWrapper.PsiBlastResult psiBlastResult);

    MutationFeatureVector createMutationFeatureVector(MutationJob mutationJob,
                                                      ChainIdentifier chainIdentifier,
                                                      ResidueIdentifier residueToMutate,
                                                      AminoAcid.Family mutationTarget);
}

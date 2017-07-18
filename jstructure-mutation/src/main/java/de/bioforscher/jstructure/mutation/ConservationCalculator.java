package de.bioforscher.jstructure.mutation;

/**
 * Extract a sequence conversation profile from a multiple-sequence alignment.
 * Created by bittrich on 7/13/17.
 */
public interface ConservationCalculator {
    /**
     * Extract a conservation profile for a given mutation job.
     * @param mutationJob the instance to process
     */
    void extractConservationProfile(MutationJob mutationJob);
}

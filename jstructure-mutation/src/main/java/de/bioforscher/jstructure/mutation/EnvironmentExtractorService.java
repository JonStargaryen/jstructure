package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;

import java.util.List;

/**
 * For a given interesting mutation, extract its surroundings and use them to find potentially similar ones. The idea is
 * that the number of potentially similar hits is exceedingly high for sequence, structure or interaction only
 * approaches. Only the combination of 3 (or more) levels will yield meaningful results, while addressing also the
 * degree of variance present in protein structures.
 * Created by bittrich on 7/6/17.
 */
public class EnvironmentExtractorService {
    private static final int WINDOW_SIZE = 3;
    private static final double INTERACTION_CUTOFF = 8.0;

    public Environment extractEnvironment(Protein protein, String chainId, int residueNumber) {
        Chain chain = protein.select()
                .chainName(chainId)
                .asChain();
        AminoAcid centralGroup = (AminoAcid) chain.select()
                .residueNumber(residueNumber)
                .asGroup();

        // evaluate sequence neighbor criteria
        GroupContainer sequentialNeighbors = chain.select()
                .residueNumber(new IntegerRange(residueNumber - WINDOW_SIZE, residueNumber + WINDOW_SIZE))
                .asGroupContainer();
        String sequence = sequentialNeighbors.getAminoAcidSequence();

        // search for spatial neighbors
        GroupContainer spatialNeighbors = chain.select()
                .aminoAcids()
                .groupDistance(centralGroup.getCa(), INTERACTION_CUTOFF)
                .asGroupContainer();

        // determine residue interaction pattern
        PLIPInteractionContainer interactionNeighbors = chain.getFeatureContainer()
                .getFeature(PLIPInteractionContainer.class)
                .getInteractionsFor(centralGroup);

        return new Environment(centralGroup,
                sequentialNeighbors,
                sequence,
                spatialNeighbors,
                interactionNeighbors);
    }

    public List<Pair<Protein, Environment>> searchForSimilarEnvironmentInDatabase(Environment environment) {
        return null;
    }
}

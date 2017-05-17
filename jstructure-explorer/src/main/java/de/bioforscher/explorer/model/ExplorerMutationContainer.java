package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;

import java.util.List;
import java.util.stream.Stream;

/**
 * The data structure of a mutated position.
 * Created by bittrich on 4/20/17.
 */
public class ExplorerMutationContainer {
    private static final double DISTANCE_CUTOFF = 8.0;
    private static final double DISTANCE_CUTOFF_SQUARED = DISTANCE_CUTOFF * DISTANCE_CUTOFF;
    private ExplorerMutationEnvironment original, mutation;

    public ExplorerMutationContainer() {

    }

    public ExplorerMutationContainer(Protein originalProtein, ExplorerChain originalExplorerChain, Protein mutatedProtein, ExplorerChain mutatedExplorerChain, int pos) {
        Group originalGroup = originalProtein.select().residueNumber(pos).asGroup();

        //TODO pattern to query here
//        String originalPlipXmlContent = PLIPRestServiceQuery.getPlipResults(originalGroup.getParentChain().getParentProtein().getName().split("_")[0], originalGroup.getParentChain().getChainId());
//        List<PLIPInteraction> originalPlipInteractions = PLIPParser.parse(originalGroup.getParentChain(), originalPlipXmlContent);
//
//        String mutatedPlipXmlContent = PLIPRestServiceQuery.getPlipResults(mutatedProtein.getName().split("_")[0], mutatedProtein.getChains().get(0).getChainId());
//        List<PLIPInteraction> mutatedPlipInteractions = PLIPParser.parse(mutatedProtein.getChains().get(0), mutatedPlipXmlContent);

        Chain originalGroups = (Chain) originalProtein.select()
                .customGroupPredicate(aminoAcid -> around(originalGroup, aminoAcid))
//                .customGroupPredicate(aminoAcid -> interacts(originalPlipInteractions, originalGroup, aminoAcid))
                .asGroupContainer();
        Chain mutatedGroups = (Chain) mutatedProtein.select()
                .customGroupPredicate(aminoAcid -> around(originalGroup, aminoAcid))
//                .customGroupPredicate(aminoAcid -> interacts(mutatedPlipInteractions, originalGroup, aminoAcid))
                .asGroupContainer();

        //TODO fetch other homologous environments, align them

        renumberSubStructures(originalGroups, mutatedGroups);

        this.original = new ExplorerMutationEnvironment(originalGroups, originalExplorerChain, originalProtein.getName());
        this.mutation = new ExplorerMutationEnvironment(mutatedGroups, mutatedExplorerChain, mutatedProtein.getName());
    }

    private void renumberSubStructures(Chain... groups) {
        Stream.of(groups)
                .forEach(chain -> {
                    for(int atomIndex = 1; atomIndex <= chain.getAtoms().size(); atomIndex++) {
                        chain.getAtoms().get(atomIndex - 1).setPdbSerial(atomIndex);
                    }
                });
    }

    private boolean interacts(List<PLIPInteraction> plipInteractions, Group referenceGroup, Group candidateGroup) {
        return plipInteractions.stream()
                .anyMatch(interaction -> (interaction.getPartner1().equals(referenceGroup) && interaction.getPartner2().equals(candidateGroup)) ||
                        (interaction.getPartner2().equals(referenceGroup) && interaction.getPartner1().equals(candidateGroup)));
    }

    /**
     * Select environment by some spatial threshold.
     * @param referenceGroup the group to look at
     * @param candidateGroup the group to evaluate
     * @return true if the candidateGroup is in spatial proximity to the reference
     */
    private boolean around(Group referenceGroup, Group candidateGroup) {
        return Combinatorics.cartesianProductOf(referenceGroup.getAtoms(), candidateGroup.getAtoms())
                .anyMatch(pair -> LinearAlgebra3D.distanceFast(pair.getLeft().getCoordinates(), pair.getRight().getCoordinates()) < DISTANCE_CUTOFF_SQUARED);
    }

    public ExplorerMutationEnvironment getOriginal() {
        return original;
    }

    public ExplorerMutationEnvironment getMutation() {
        return mutation;
    }
}

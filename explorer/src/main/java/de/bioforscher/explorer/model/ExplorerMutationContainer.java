package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;

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

        Chain originalGroups = (Chain) originalProtein.select()
                .customGroupPredicate(aminoAcid -> around(originalGroup, aminoAcid))
                .asGroupContainer();
        Chain mutatedGroups = (Chain) mutatedProtein.select()
                .customGroupPredicate(aminoAcid -> around(originalGroup, aminoAcid))
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

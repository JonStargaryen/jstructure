package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public interface Interaction {
    /**
     * @return true for intra-chain interactions (amino acid - amino acid)
     */
    boolean isPolymerInteraction();

    /**
     * @return true for ligand interactions (amino acid - small-molecule)
     */
    boolean isLigandInteraction();

    /**
     * @return true for directed interactions such as hydrogen bonds
     */
    boolean isDirectedInteraction();

    default Stream<Atom> interactingAtoms1() {
        return getInteractingAtoms1().stream();
    }

    List<Atom> getInteractingAtoms1();

    default Stream<Atom> interactingAtoms2() {
        return getInteractingAtoms2().stream();
    }

    List<Atom> getInteractingAtoms2();

    default Stream<Atom> allInteractingAtoms() {
        return Stream.concat(interactingAtoms1(), interactingAtoms2());
    }

    default List<Atom> getAllInteractingAtoms() {
        return allInteractingAtoms()
                .collect(Collectors.toList());
    }

    boolean containsGroup(Group group);
}

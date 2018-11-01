package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;

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

    Group getInteractingGroup1();

    Group getInteractingGroup2();

    Group getOpposingGroup(Group group);

    default double[] getCentroid() {
        return LinearAlgebra.on(interactingAtoms1().collect(StructureCollectors.toCentroid()))
                .add(interactingAtoms2().collect(StructureCollectors.toCentroid()))
                .divide(2)
                .getValue();
//        Atom  atom = Atom.builder(mapInteractionToElement(this), centroid).build();
//        if(isLigandInteraction()) {
//            Group opposingGroup = getInteractingGroup1().isLigand() ? getInteractingGroup2() : getInteractingGroup1();
//            opposingGroup.addAtom(atom);
//        }
//        return atom;
    }

//    static Element mapInteractionToElement(Interaction interaction) {
//        if(interaction instanceof HalogenBond) {
//            return Element.H;
//        } else if(interaction instanceof HydrogenBond) {
//            return Element.O;
//        } else if(interaction instanceof HydrophobicInteraction) {
//            return Element.C;
//        } else if(interaction instanceof MetalComplex) {
//            return Element.Fe;
//        } else if(interaction instanceof PiCationInteraction) {
//            return Element.Ca;
//        } else if(interaction instanceof PiStackingInteraction) {
//            return Element.K;
//        } else if(interaction instanceof SaltBridge) {
//            return Element.Cl;
//        } else {
//            return Element.Ar;
//        }
//    }
}

package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.mathematics.CoordinateManipulations;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The abstract capabilities of algorithms which merge several observations of fragments or motifs into 1 consensus
 * representation.
 * Created by S on 05.12.2016.
 */
public abstract class AbstractConsensusComposer {
    static final Logger logger = LoggerFactory.getLogger(AbstractConsensusComposer.class);

    private static Atom mergeAtomPair(Pair<Atom, Atom> atomPair) {
        Atom atom = atomPair.getLeft();
        double[] left = atom.getCoordinates();
        double[] right = atomPair.getRight().getCoordinates();
        atom.setCoordinates(LinearAlgebra3D.divide(LinearAlgebra3D.add(left, right), 2));
        return atom;
    }

    static AtomContainer mergeContainerPair(AtomContainer reference, AtomContainer candidate) {
        // get intersecting pairs of atoms wrapped in an atom container
        Pair<AtomContainer, AtomContainer> containerPair =
                CoordinateManipulations.comparableAtomContainerPair(reference, candidate);

        return Combinatorics.sequentialPairsOf(containerPair.getLeft().getAtoms(), containerPair.getRight().getAtoms())
                .map(AbstractConsensusComposer::mergeAtomPair)
                .collect(StructureCollectors.toAtomContainer());
    }
}

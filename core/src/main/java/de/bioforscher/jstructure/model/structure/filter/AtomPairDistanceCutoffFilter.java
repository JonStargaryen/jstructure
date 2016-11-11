package de.bioforscher.jstructure.model.structure.filter;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;

import java.util.function.Predicate;

/**
 * A filter which returns <code>true</code> when the euclidean distance between a {@link Pair} of atoms is smaller
 * than a given threshold.
 */
public class AtomPairDistanceCutoffFilter implements Predicate<Pair<Atom, Atom>> {
    private final double squaredDistanceCutoff;

    public AtomPairDistanceCutoffFilter(double distanceCutoff) {
        // square for faster computations
        this.squaredDistanceCutoff = distanceCutoff * distanceCutoff;
    }

    @Override
    public boolean test(Pair<Atom, Atom> atomPair) {
        return LinearAlgebra3D.distanceFast(atomPair.getLeft().getCoordinates(),
                atomPair.getRight().getCoordinates()) < squaredDistanceCutoff;
    }
}

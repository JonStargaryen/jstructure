package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;

import java.util.function.Predicate;

/**
 * A select which returns <code>true</code> for each atom of a collection whose euclidean distance to a reference
 * atom is smaller than a given threshold.
 */
@Deprecated
public class AtomDistanceCutoffFilter implements Predicate<Atom> {
    private final Atom referenceAtom;
    private final double[] referenceCoordinates;
    private final double squaredDistanceCutoff;

    public AtomDistanceCutoffFilter(Atom referenceAtom, double distanceCutoff) {
        this.referenceAtom = referenceAtom;
        this.referenceCoordinates = this.referenceAtom.getCoordinates();
        // square for faster computations
        this.squaredDistanceCutoff = distanceCutoff * distanceCutoff;
    }

    @Override
    public boolean test(Atom atomToTest) {
        // return true when distance is smaller than threshold, but ensure atoms are not equal
        return !atomToTest.equals(referenceAtom) && LinearAlgebra3D.distanceFast(referenceCoordinates,
                atomToTest.getCoordinates()) < squaredDistanceCutoff;
    }
}

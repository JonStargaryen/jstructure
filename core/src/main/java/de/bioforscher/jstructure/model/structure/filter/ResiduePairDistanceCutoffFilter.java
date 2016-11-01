package de.bioforscher.jstructure.model.structure.filter;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;

import java.util.function.Predicate;

public class ResiduePairDistanceCutoffFilter implements Predicate<Pair<Residue, Residue>> {
    private final double squaredDistanceCutoff;
    private final RepresentationScheme representationScheme;

    public ResiduePairDistanceCutoffFilter(double distanceCutoff, RepresentationScheme representationScheme) {
        this.squaredDistanceCutoff = distanceCutoff * distanceCutoff;
        this.representationScheme = representationScheme;
    }

    @Override
    public boolean test(Pair<Residue, Residue> residuePair) {
        return LinearAlgebra3D.distanceFast(representationScheme.determineRepresentingAtom(residuePair.getFirst()).getCoordinates(),
                representationScheme.determineRepresentingAtom(residuePair.getSecond()).getCoordinates()) < squaredDistanceCutoff;
    }
}
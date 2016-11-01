package de.bioforscher.jstructure.model.structure.filter;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.scheme.CentroidRepresentationScheme;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;

import java.util.function.Predicate;

public class GroupPairDistanceCutoffFilter implements Predicate<Pair<Group, Group>> {
    private final double squaredDistanceCutoff;
    private final RepresentationScheme representationScheme;

    public GroupPairDistanceCutoffFilter(double distanceCutoff) {
        this.squaredDistanceCutoff = distanceCutoff * distanceCutoff;
        this.representationScheme = new CentroidRepresentationScheme();
    }

    @Override
    public boolean test(Pair<Group, Group> groupPair) {
        return LinearAlgebra3D.distanceFast(representationScheme.determineRepresentingAtom(groupPair.getFirst()).getCoordinates(),
                representationScheme.determineRepresentingAtom(groupPair.getSecond()).getCoordinates()) <
                squaredDistanceCutoff;
    }
}
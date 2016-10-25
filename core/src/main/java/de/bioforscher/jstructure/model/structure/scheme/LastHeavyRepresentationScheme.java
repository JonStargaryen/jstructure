package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Residue;

/**
 * Represents a residue by the non-hydrogen atom most distant to the alpha carbon.
 * Created by S on 05.10.2016.
 */
public class LastHeavyRepresentationScheme implements RepresentationScheme {
    private static final RepresentationScheme alphaCarbonRepresentationScheme = new AlphaCarbonRepresentationScheme();

    @Override
    public Atom determineRepresentingAtom(Group group) {
        if(!(group instanceof Residue)) {
            throw new IllegalArgumentException(getClass().getSimpleName() +
                    " can only be employed to residue objects use CentroidScheme instead");
        }
        final Atom alphaCarbon = alphaCarbonRepresentationScheme.determineRepresentingAtom(group);
        return group.nonHydrogenAtoms()
                    .max((a1, a2) -> LinearAlgebra3D.distanceFast(alphaCarbon.getCoordinates(), a1.getCoordinates()) <
                        LinearAlgebra3D.distanceFast(alphaCarbon.getCoordinates(), a2.getCoordinates()) ? 1 : -1).get();
    }
}

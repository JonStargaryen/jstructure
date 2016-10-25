package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Residue;

/**
 * Represents a group by its alpha carbon atom. The centroid of all atoms is used as fallback.
 * Created by S on 02.06.2016.
 */
public class AlphaCarbonRepresentationScheme implements RepresentationScheme {
    @Override
    public Atom determineRepresentingAtom(Group group) {
        if(!(group instanceof Residue)) {
            throw new IllegalArgumentException(getClass().getSimpleName() +
                    " can only be employed to residue objects use CentroidScheme instead");
        }
        return ((Residue) group).alphaCarbon().orElse(new Atom(CoordinateUtils.centroid(group.atoms())));
    }
}

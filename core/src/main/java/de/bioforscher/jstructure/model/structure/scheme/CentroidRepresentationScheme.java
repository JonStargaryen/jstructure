package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

/**
 * Represents a group of atoms by the centroid of all atoms.
 * Created by S on 05.10.2016.
 */
public class CentroidRepresentationScheme implements RepresentationScheme {
    @Override
    public Atom determineRepresentingAtom(Group group) {
        return new Atom(CoordinateUtils.centroid(group));
    }
}

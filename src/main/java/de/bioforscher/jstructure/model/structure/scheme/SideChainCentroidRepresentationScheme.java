package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Residue;

/**
 * Created by S on 05.10.2016.
 */
public class SideChainCentroidRepresentationScheme implements RepresentationScheme {
    @Override
    public Atom determineRepresentingAtom(Group group) {
        if(!(group instanceof Residue)) {
            throw new IllegalArgumentException(getClass().getSimpleName() + " can only be employed to residue objects use CentroidScheme instead");
        }
        return new Atom(CoordinateUtils.centroid(((Residue) group).sideChainAtoms()));
    }
}

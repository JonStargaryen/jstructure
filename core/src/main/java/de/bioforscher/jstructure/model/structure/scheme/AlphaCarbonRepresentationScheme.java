package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.CoordinateManipulations;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.selection.Selection;

/**
 * Represents a group by its alpha carbon atom. The centroid of all atoms is used as fallback.
 * Created by S on 02.06.2016.
 */
public class AlphaCarbonRepresentationScheme implements RepresentationScheme {
    @Override
    public Atom determineRepresentingAtom(Group group) {
        if(!group.isAminoAcid()) {
            throw new IllegalArgumentException(getClass().getSimpleName() +
                    " can only be employed to getResidue objects use CentroidScheme instead");
        }
        return Selection.on(group)
                .alphaCarbonAtoms()
                .asOptionalAtom()
                .orElse(new Atom(CoordinateManipulations.centroid(group)));
    }
}

package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.CoordinateManipulations;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.selection.Selection;

/**
 * Represents an amino acid by its beta carbon atom. The centroid of all atoms will be used as fallback (i.e. primarily
 * the case for glycine).
 * Created by S on 02.06.2016.
 */
public class BetaCarbonRepresentationScheme implements RepresentationScheme {
    @Override
    public Atom determineRepresentingAtom(Group group) {
        if(!group.isAminoAcid()) {
            throw new IllegalArgumentException(getClass().getSimpleName() +
                    " can only be employed to amino acids - use CentroidScheme instead");
        }
        return Selection.on(group)
                .betaCarbonAtoms()
                .asOptionalAtom()
                .orElse(new Atom(CoordinateManipulations.centroid(group)));
    }
}

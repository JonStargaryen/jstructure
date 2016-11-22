package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

/**
 * Created by S on 08.11.2016.
 */
public class SideChainCentroidRepresentationScheme implements RepresentationScheme {
    @Override
    public Atom determineRepresentingAtom(Group group) {
        if(!group.isAminoAcid()) {
            throw new IllegalArgumentException(getClass().getSimpleName() +
                    " can only be employed to residue objects use CentroidScheme instead");
        }
        AtomContainer sideChainAtoms = Selection.on(group)
                .sideChainAtoms()
                .asAtomContainer();
        return new Atom(CoordinateUtils.centroid(sideChainAtoms));
    }
}

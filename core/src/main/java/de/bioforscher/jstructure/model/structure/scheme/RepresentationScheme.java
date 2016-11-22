package de.bioforscher.jstructure.model.structure.scheme;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

/**
 * Defines ways to represent residues by a single atom. This can be used e.g. for distance calculations based on alpha
 * carbons only.
 * Created by S on 05.10.2016.
 */
public interface RepresentationScheme {
    /**
     * determine the coordinates best fit to represent this getResidue - this is a key step in calculating distances
     * between residues pair and such
     * @param group the group to be processed
     * @return an atom (essentially we only are interested in its coordinates) representing the given {@link Group}
     * by the rules defined by the impl
     */
    Atom determineRepresentingAtom(Group group);
}

package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import de.bioforscher.jstructure.model.structure.filter.AtomPairDistanceCutoffFilter;

import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Defines the capability to convert itself to <tt>ATOM</tt> records according to <tt>PDB</tt> format.
 * Created by S on 28.09.2016.
 */
public interface AtomContainer extends Container {
    /**
     * Access to all atoms of this container.
     * @return a stream of all atoms in this container
     */
    Stream<Atom> atoms();

    /**
     * Discards all atoms from this container.
     */
    void clear();

    /**
     * @return a stream of all non-hydrogen atoms
     * @see AtomNameFilter#HYDROGEN_FILTER
     */
    default Stream<Atom> nonHydrogenAtoms() {
        return atoms().filter(AtomNameFilter.HYDROGEN_FILTER.negate());
    }

    /**
     * Returns all atom pairs which are closer to each other than defined by some distance threshold.
     * @param distanceCutoff a distance threshold in Angstrom
     * @return a stream of all atom pairs in contact
     */
    default Stream<Pair<Atom, Atom>> atomPairsInContact(double distanceCutoff) {
        final AtomPairDistanceCutoffFilter distanceCutoffFilter = new AtomPairDistanceCutoffFilter(distanceCutoff);
        return Pair.uniquePairsOf(atoms().collect(Collectors.toList())).filter(distanceCutoffFilter);
    }

    /**
     * Returns all atom pairs.
     * @return a stream of all atom pairs
     */
    default Stream<Pair<Atom, Atom>> atomPairs() {
        return Pair.uniquePairsOf(atoms().collect(Collectors.toList()));
    }

    /**
     * Deletes all pseudo atoms originating by internal computations. They are identified by a pdbSerial of
     * {@link Integer#MIN_VALUE}.
     */
    void clearPseudoAtoms();
}
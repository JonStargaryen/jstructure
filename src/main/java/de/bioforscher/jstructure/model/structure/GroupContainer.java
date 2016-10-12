package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import de.bioforscher.jstructure.model.structure.scheme.AlphaCarbonRepresentationScheme;

import java.util.stream.Stream;

/**
 * Specifies the capabilities of a group container (mostly a {@link Chain}).
 * Created by S on 30.09.2016.
 */
public interface GroupContainer extends AtomContainer {
    /**
     * Access to all groups associated to this container.
     * @return a stream of groups
     */
    Stream<Group> groups();

    Stream<Residue> residues();

    /**
     * Returns all residue pairs which are closer to each other than defined by some distance threshold.
     * @param distanceCutoff a distance threshold in Angstrom
     * @return a stream of all residue pairs in contact
     * @see CoordinateUtils#atomPairsInContact(Stream, double)
     */
    default Stream<Pair<Residue>> residuesInContact(double distanceCutoff) {
        return CoordinateUtils.residuePairsInContact(residues(), distanceCutoff, new AlphaCarbonRepresentationScheme());
    }

    /**
     * @return a stream of non-hydrogen atoms
     * @see AtomNameFilter#HYDROGEN_FILTER
     */
    default Stream<Atom> nonHydrogenAtoms() {
        return atoms().filter(AtomNameFilter.HYDROGEN_FILTER.negate());
    }

    /**
     * @return a stream of all backbone hydrogen atoms
     * @see AtomNameFilter#HYDROGEN_FILTER
     */
    default Stream<Atom> backboneHydrogens() {
        return atoms().filter(AtomNameFilter.HYDROGEN_FILTER);
    }

    default Stream<Atom> backboneCarbon() {
        return atoms().filter(AtomNameFilter.C_ATOM_FILTER);
    }

    default Stream<Atom> backboneOxygen() {
        return atoms().filter(AtomNameFilter.O_ATOM_FILTER);
    }

    default Stream<Atom> backboneNitrogen() {
        return atoms().filter(AtomNameFilter.N_ATOM_FILTER);
    }

    /**
     * @return a stream of alpha carbons
     * @see AtomNameFilter#CA_ATOM_FILTER
     */
    default Stream<Atom> alphaCarbons() {
        return atoms().filter(AtomNameFilter.CA_ATOM_FILTER);
    }

    default Stream<Atom> betaCarbons() {
        return atoms().filter(AtomNameFilter.CB_ATOM_FILTER);
    }
}
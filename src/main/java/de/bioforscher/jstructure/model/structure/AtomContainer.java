package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.Pair;

import java.util.List;
import java.util.function.Predicate;
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
     * @return a stream of backbone atoms
     * @see AtomNameFilter#BACKBONE_ATOM_FILTER
     */
    default Stream<Atom> backboneAtoms() {
        return atoms().filter(AtomNameFilter.BACKBONE_ATOM_FILTER);
    }

    /**
     * @return a stream of alpha carbons
     * @see AtomNameFilter#CA_ATOM_FILTER
     */
    default Stream<Atom> alphaCarbons() {
        return atoms().filter(AtomNameFilter.CA_ATOM_FILTER);
    }

    /**
     * @return a stream of non-hydrogen atoms
     * @see AtomNameFilter#NON_HYDROGEN_FILTER
     */
    default Stream<Atom> nonHydrogenAtoms() {
        return atoms().filter(AtomNameFilter.NON_HYDROGEN_FILTER);
    }

    /**
     * @return a stream of side chain atoms
     * @see AtomNameFilter#SIDE_CHAIN_ATOM_FILTER
     */
    default Stream<Atom> sideChainAtoms() {
        return atoms().filter(AtomNameFilter.SIDE_CHAIN_ATOM_FILTER);
    }

    /**
     *
     * @param distanceCutoff a distance threshold in Angstrom
     * @return a stream of all atom pairs in contact
     * @see CoordinateUtils#atomPairsInContact(Stream, double)
     */
    default Stream<Pair<Atom>> atomPairsInContact(double distanceCutoff) {
        return CoordinateUtils.atomPairsInContact(atoms(), distanceCutoff);
    }

    /**
     * Tests whether an atom is of a certain <tt>PDB</tt> name. This allows for filtering for CA or backbone atoms.
     */
    class AtomNameFilter implements Predicate<Atom> {
        private List<String> acceptedAtomNames;
        private boolean negate;
        /**
         * Retains all backbone atoms (i.e. 'N', 'CA', 'C' or 'O').
         */
        public final static AtomNameFilter BACKBONE_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.BACKBONE_ATOM_NAMES, false);
        /**
         * Retains all atoms which are side chain atoms (i.e. part of the atom names occurring for the 20 standard amino
         * acid, but neither 'N', 'CA', 'C' or 'O'.
         */
        public final static AtomNameFilter SIDE_CHAIN_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.SIDECHAIN_ATOM_NAMES, false);
        /**
         * Retains all atoms whose {@link Atom#getName()} equals the alpha carbon name 'CA'.
         */
        public final static AtomNameFilter CA_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.CA_ATOM_NAMES, false);
        /**
         * This filter retains all non-hydrogen atoms (i.e. their element does not equal {@link Element#H},
         * {@link Element#D} or {@link Element#T}.
         */
        public final static AtomNameFilter NON_HYDROGEN_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.H_ATOM_NAME, true);
        /**
         * This filter retains all atom names.
         */
        public final static AtomNameFilter ALL_ATOM_FILTER = new AtomNameFilter(null, true) {
            @Override
            public boolean test(Atom atom) {
                return true;
            }
        };

        public AtomNameFilter(final List<String> acceptedAtomNames, boolean negate) {
            this.acceptedAtomNames = acceptedAtomNames;
            this.negate = negate;
        }

        @Override
        public boolean test(Atom atom) {
            return negate != this.acceptedAtomNames.contains(atom.getName());
        }
    }
}

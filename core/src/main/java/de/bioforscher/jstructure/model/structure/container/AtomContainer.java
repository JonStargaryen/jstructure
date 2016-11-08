package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import de.bioforscher.jstructure.model.structure.filter.AtomPairDistanceCutoffFilter;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Defines the capability to convert itself to <tt>ATOM</tt> records according to <tt>PDB</tt> format.
 * Created by S on 28.09.2016.
 */
public interface AtomContainer extends Container {
    List<Atom> getAtoms();

    /**
     * Access to all atoms of this container.
     * @return a stream of all atoms in this container
     */
    default Stream<Atom> atoms() {
        return getAtoms().stream();
    }

    /**
     * Discards all atoms from this container.
     */
    default void clear() {
        getAtoms().clear();
    }

    /**
     * Deletes all pseudo atoms originating by internal computations. They are identified by a pdbSerial of
     * {@link Integer#MIN_VALUE}.
     */
    default void clearPseudoAtoms() {
        getAtoms().removeIf(Atom::isVirtual);
    }

    /**
     * @return a stream of all non-hydrogen atoms
     * @see AtomNameFilter#HYDROGEN_FILTER
     */
    default Stream<Atom> nonHydrogenAtoms() {
        return atoms().filter(AtomNameFilter.HYDROGEN_FILTER.negate());
    }

    default String composePDBRecord() {
        return atoms().map(Atom::composePDBRecord)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    /**
     * Returns all atom pairs which are closer to each other than defined by some distance threshold.
     * @param distanceCutoff a distance threshold in Angstrom
     * @return a stream of all atom pairs in contact
     */
    default Stream<Pair<Atom, Atom>> atomPairsInContact(double distanceCutoff) {
        return atomPairs().filter(new AtomPairDistanceCutoffFilter(distanceCutoff));
    }

    /**
     * Returns all atom pairs.
     * @return a stream of all atom pairs
     */
    default Stream<Pair<Atom, Atom>> atomPairs() {
        return Pair.uniquePairsOf(getAtoms());
    }

    static AtomContainer of(Stream<Atom> atomStream) {
        return atomStream.collect(AtomContainerConsumer::new,
                AtomContainerConsumer::accept,
                AtomContainerConsumer::combine)
            .getContainer();
    }

    class AtomContainerConsumer implements Consumer<Atom> {
        List<Atom> atoms;

        AtomContainerConsumer() {
            this.atoms = new ArrayList<>();
        }

        AtomContainer getContainer() {
            return new BasicAtomContainer(atoms);
        }

        @Override
        public void accept(Atom atom) {
            atoms.add(atom);
        }

        void combine(AtomContainerConsumer other) {
            atoms.addAll(other.atoms);
        }
    }

    class BasicAtomContainer implements AtomContainer {
        private List<Atom> atoms;

        BasicAtomContainer(List<Atom> atomList) {
            this.atoms = atomList;
        }

        @Override
        public List<Atom> getAtoms() {
            return atoms;
        }

        @Override
        public Map<Enum, Object> getFeatureMap() {
            throw new UnsupportedOperationException();
        }
    }
}
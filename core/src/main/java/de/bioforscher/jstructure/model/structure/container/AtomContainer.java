package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Defines the capability to convert itself to <tt>ATOM</tt> records according to <tt>PDB</tt> format.
 * Created by S on 28.09.2016.
 */
public interface AtomContainer extends StructureContainer {
    List<Atom> getAtoms();

    /**
     * Access to all atoms of this container.
     * @return a select of all atoms in this container
     */
    default Stream<Atom> atoms() {
        return getAtoms().stream();
    }

    default String composePDBRecord() {
        return atoms().map(Atom::composePDBRecord)
                      .collect(Collectors.joining(System.lineSeparator()));
    }

    /**
     * Discards all atoms from this container.
     */
    default void clearAtoms() {
        getAtoms().clear();
    }

    /**
     * Deletes all pseudo atoms originating by internal computations.
     */
    default void clearPseudoAtoms() {
        getAtoms().removeIf(Atom::isVirtual);
    }

    static AtomContainer of(Atom... atoms) {
        return new Group(Arrays.asList(atoms));
    }
}
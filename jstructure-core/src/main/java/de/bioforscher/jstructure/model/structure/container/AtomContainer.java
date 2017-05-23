package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.Algebrable;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.selection.Selectable;

import java.lang.reflect.InvocationTargetException;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Defines the capability to convert itself to <tt>ATOM</tt> records according to <tt>PDB</tt> format.
 * Created by S on 28.09.2016.
 */
public interface AtomContainer extends StructureContainer, Selectable, Algebrable<LinearAlgebra.AtomContainerLinearAlgebra> {
    /**
     * Never manipulate the returned collection as it is not guaranteed the actually modify the internal list(s).
     * @return all associated atoms
     */
    List<Atom> getAtoms();

    /**
     * Access to all atoms of this container.
     * @return a select of all atoms in this container
     */
    default Stream<Atom> atoms() {
        return getAtoms().stream();
    }

    default String getPdbRepresentation() {
        // TODO the current approach ignoring TER records
        return atoms()
                .sorted(Comparator.comparingInt(Atom::getPdbSerial))
                .map(Atom::getPdbRepresentation)
                .collect(Collectors.joining(System.lineSeparator(),
                        "",
                        System.lineSeparator()));
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

    default AtomContainer createCopy() {
        try {
            return getClass().getConstructor(getClass()).newInstance(this);
        } catch (NoSuchMethodException | IllegalAccessException | InstantiationException | InvocationTargetException e) {
            throw new UnsupportedOperationException(e);
        }
    }
}
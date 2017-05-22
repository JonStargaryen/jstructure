package de.bioforscher.jstructure.model.structure;

/**
 * Specifies the ability to convert itself to an <tt>ATOM</tt> record representation in <tt>PDB</tt> format.
 * Created by S on 30.09.2016.
 */
public interface AtomRecordWriter {
    /**
     * Composes the <tt>PDB</tt> representation of this element by creating an <tt>ATOM</tt> record for each
     * {@link Atom} within this container.
     * @return an <tt>ATOM</tt> for each {@link Atom} associated to this element, separated by
     * {@link System#lineSeparator()}
     */
    String getPdbRepresentation();
}
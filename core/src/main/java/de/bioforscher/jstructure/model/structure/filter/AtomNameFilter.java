package de.bioforscher.jstructure.model.structure.filter;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;

import java.util.List;
import java.util.Objects;
import java.util.function.Predicate;

/**
 * Tests whether an atom is of a certain <tt>PDB</tt> name. This allows for filtering for CA or backbone atoms.
 */
public class AtomNameFilter implements Predicate<Atom> {
    private final List<String> acceptedAtomNames;
    /**
     * Retains all backbone atoms (i.e. 'N', 'CA', 'C' or 'O').
     */
    public final static AtomNameFilter BACKBONE_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.BACKBONE_ATOM_NAMES);
    /**
     * Retains all atoms which are side chain atoms (i.e. part of the atom names occurring for the 20 standard amino
     * acid, but neither 'N', 'CA', 'C' or 'O'.
     */
    public final static AtomNameFilter SIDE_CHAIN_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.SIDECHAIN_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the alpha carbon name 'CA'.
     */
    public final static AtomNameFilter CA_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.CA_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the backbone carbon name 'C'.
     */
    public final static AtomNameFilter C_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.C_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the backbone nitrogen name 'N'.
     */
    public final static AtomNameFilter N_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.N_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the backbone oxygen name 'O'.
     */
    public final static AtomNameFilter O_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.O_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the beta carbon name 'CB'.
     */
    public final static AtomNameFilter CB_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.CA_ATOM_NAMES);

    /**
     * This filter retains all hydrogen atoms (i.e. their element does equal {@link Element#H},
     * {@link Element#D} or {@link Element#T}.
     */
    public final static AtomNameFilter HYDROGEN_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.H_ATOM_NAME);

    public AtomNameFilter(final List<String> acceptedAtomNames) {
        this.acceptedAtomNames = Objects.requireNonNull(acceptedAtomNames);
    }

    @Override
    public boolean test(Atom atom) {
        return acceptedAtomNames.contains(atom.getName());
    }
}

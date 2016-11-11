package de.bioforscher.jstructure.model.structure.filter;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Tests whether an atom is of a certain <tt>PDB</tt> name. This allows for filtering for CA or backbone atoms.
 */
public class AtomNameFilter implements Predicate<Atom> {
    /**
     * Retains all backbone atoms (i.e. 'N', 'CA', 'C' or 'O').
     */
    public static final AtomNameFilter BACKBONE_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.BACKBONE_ATOM_NAMES);
    /**
     * Retains all atoms which are side getChain atoms (i.e. part of the atom names occurring for the 20 standard amino
     * acid, but neither 'N', 'CA', 'C' or 'O'.
     */
    public static final AtomNameFilter SIDE_CHAIN_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.SIDECHAIN_ATOM_NAMES);
    /**
     * Retains all atoms.
     */
    public static final AtomNameFilter ALL_ATOM_FILTER = new AtomNameFilter(Stream.concat(AminoAcid.ATOM_NAMES.BACKBONE_ATOM_NAMES.stream(),
            AminoAcid.ATOM_NAMES.SIDECHAIN_ATOM_NAMES.stream()).collect(Collectors.toList()));
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the alpha carbon name 'CA'.
     */
    public static final AtomNameFilter CA_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.CA_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the backbone carbon name 'C'.
     */
    public static final AtomNameFilter C_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.C_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the backbone nitrogen name 'N'.
     */
    public static final AtomNameFilter N_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.N_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the backbone oxygen name 'O'.
     */
    public static final AtomNameFilter O_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.O_ATOM_NAMES);
    /**
     * Retains all atoms whose {@link Atom#getName()} equals the beta carbon name 'CB'.
     */
    public static final AtomNameFilter CB_ATOM_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.CA_ATOM_NAMES);

    /**
     * This filter retains all hydrogen atoms (i.e. their element does equal {@link Element#H},
     * {@link Element#D} or {@link Element#T}.
     */
    public static final AtomNameFilter HYDROGEN_FILTER = new AtomNameFilter(AminoAcid.ATOM_NAMES.H_ATOM_NAME);

    public static final AtomNameFilter O3PRIME_FILTER = new AtomNameFilter(Collections.singletonList("O3'"));

    public static final AtomNameFilter O5PRIME_FILTER = new AtomNameFilter(Collections.singletonList("O5'"));

    public static final AtomNameFilter PHOSPHATE_FILTER = new AtomNameFilter(Collections.singletonList("P"));

    private final List<String> acceptedAtomNames;

    public List<String> getAcceptedAtomNames() {
        return acceptedAtomNames;
    }

    public AtomNameFilter(final List<String> acceptedAtomNames) {
        this.acceptedAtomNames = Objects.requireNonNull(acceptedAtomNames);
    }

    @Override
    public boolean test(Atom atom) {
        return acceptedAtomNames.contains(atom.getName());
    }
}

package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;

import java.util.Optional;
import java.util.stream.Stream;

/**
 * Represents one amino acid within a {@link Protein} and is identified by a getResidue number, e.g. ALA-103. Is composed
 * of {@link Atom} objects.
 * Created by S on 27.09.2016.
 */
public class Residue extends Group {
    private AminoAcid aminoAcid;

    /**
     * The constructor of residues.
     * @param pdbName the amino acid this getResidue represents as String
     * @param residueNumber a {@link Chain}-wide unique number to identify this getResidue
     * @see AminoAcid#valueOfIgnoreCase(String)
     */
    public Residue(String pdbName, int residueNumber) {
        super(pdbName, residueNumber);
        this.aminoAcid = AminoAcid.valueOfIgnoreCase(pdbName);
    }

    /**
     * Returns which of the 20 canonical amino acids this getResidue represents.
     * @return the type of this getResidue
     * @see AminoAcid
     */
    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }

    /**
     * Sets the amino acid for this getResidue.
     * @param aminoAcid the new value
     */
    public void setAminoAcid(AminoAcid aminoAcid) {
        this.aminoAcid = aminoAcid;
    }

    public Atom getAlphaCarbon() {
        return getAtomByName(AtomNameFilter.CA_ATOM_FILTER);
    }

    public Optional<Atom> findAlphaCarbon() {
        return tryToGetAtomByName(AtomNameFilter.CA_ATOM_FILTER);
    }

    public Atom getBackboneCarbon() {
        return getAtomByName(AtomNameFilter.C_ATOM_FILTER);
    }

    public Optional<Atom> findBackboneCarbon() {
        return tryToGetAtomByName(AtomNameFilter.C_ATOM_FILTER);
    }

    public Atom getBackboneOxygen() {
        return getAtomByName(AtomNameFilter.O_ATOM_FILTER);
    }

    public Optional<Atom> findBackboneOxygen() {
        return tryToGetAtomByName(AtomNameFilter.O_ATOM_FILTER);
    }

    public Atom getBetaCarbon() {
        return getAtomByName(AtomNameFilter.CB_ATOM_FILTER);
    }

    public Optional<Atom> findBetaCarbon() {
        return tryToGetAtomByName(AtomNameFilter.CB_ATOM_FILTER);
    }

    public Atom getBackboneNitrogen() {
        return getAtomByName(AtomNameFilter.N_ATOM_FILTER);
    }

    public Optional<Atom> findBackboneNitrogen() {
        return tryToGetAtomByName(AtomNameFilter.N_ATOM_FILTER);
    }

    public Atom getBackboneHydrogen() { return getAtomByName(AtomNameFilter.HYDROGEN_FILTER); }

    public Optional<Atom> findBackboneHydrogen() {
        return tryToGetAtomByName(AtomNameFilter.HYDROGEN_FILTER);
    }

    /**
     * @return a stream of backbone atoms
     * @see AtomNameFilter#BACKBONE_ATOM_FILTER
     */
    public Stream<Atom> backboneAtoms() {
        return atoms().filter(AtomNameFilter.BACKBONE_ATOM_FILTER);
    }

    /**
     * @return a stream of side getChain atoms
     * @see AtomNameFilter#SIDE_CHAIN_ATOM_FILTER
     */
    public Stream<Atom> sideChainAtoms() {
        return atoms().filter(AtomNameFilter.SIDE_CHAIN_ATOM_FILTER);
    }
}
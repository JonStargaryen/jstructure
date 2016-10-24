package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;

import java.util.Optional;
import java.util.stream.Stream;

/**
 * Represents one amino acid within a {@link Protein} and is identified by a residue number, e.g. ALA-103. Is composed
 * of {@link Atom} objects.
 * Created by S on 27.09.2016.
 */
public class Residue extends Group {
    private AminoAcid aminoAcid;

    /**
     * The constructor of residues.
     * @param aminoAcid the amino acid this residue represents as Object
     * @param residueNumber a {@link Chain}-wide unique number to identify this residue
     * @see AminoAcid
     */
    public Residue(AminoAcid aminoAcid, int residueNumber) {
        super(aminoAcid.getThreeLetterCode(), residueNumber);
        this.aminoAcid = aminoAcid;
    }

    /**
     * The constructor of residues.
     * @param aminoAcidName the amino acid this residue represents as String
     * @param residueNumber a {@link Chain}-wide unique number to identify this residue
     * @see AminoAcid#valueOfIgnoreCase(String)
     */
    public Residue(String aminoAcidName, int residueNumber) {
        this(AminoAcid.valueOfIgnoreCase(aminoAcidName), residueNumber);
    }

    /**
     * Returns which of the 20 canonical amino acids this residue represents.
     * @return the type of this residue
     * @see AminoAcid
     */
    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " name='" + aminoAcid + "' resNum='" + residueNumber + "' size='" + atoms.size() + "'";
    }

    public Optional<Atom> alphaCarbon() {
        return atoms().filter(AtomNameFilter.CA_ATOM_FILTER).findFirst();
    }

    public Optional<Atom> backboneCarbon() {
        return atoms().filter(AtomNameFilter.C_ATOM_FILTER).findFirst();
    }

    public Optional<Atom> backboneOxygen() {
        return atoms().filter(AtomNameFilter.O_ATOM_FILTER).findFirst();
    }

    public Optional<Atom> betaCarbon() {
        return atoms().filter(AtomNameFilter.CB_ATOM_FILTER).findFirst();
    }

    public Optional<Atom> backboneNitrogen() {
        return atoms().filter(AtomNameFilter.N_ATOM_FILTER).findFirst();
    }

    public Optional<Atom> backboneHydrogen() { return atoms().filter(AtomNameFilter.BACKBONE_ATOM_FILTER)
            .filter(AtomNameFilter.HYDROGEN_FILTER).findFirst(); }


    /**
     * @return a stream of backbone atoms
     * @see AtomNameFilter#BACKBONE_ATOM_FILTER
     */
    public Stream<Atom> backboneAtoms() {
        return atoms().filter(AtomNameFilter.BACKBONE_ATOM_FILTER);
    }

    /**
     * @return a stream of side chain atoms
     * @see AtomNameFilter#SIDE_CHAIN_ATOM_FILTER
     */
    public Stream<Atom> sideChainAtoms() {
        return atoms().filter(AtomNameFilter.SIDE_CHAIN_ATOM_FILTER);
    }
}
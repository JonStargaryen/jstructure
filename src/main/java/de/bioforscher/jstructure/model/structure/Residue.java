package de.bioforscher.jstructure.model.structure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents one amino acid within a {@link Protein} and is identified by a residue number, e.g. ALA-103. Is composed
 * of {@link Atom} objects.
 * Created by S on 27.09.2016.
 */
public class Residue implements AtomContainer, AtomRecordWriter {
    private List<Atom> atoms;
    private int residueNumber;
    private AminoAcid aminoAcid;
    /**
     * Handle to the container element.
     */
    private Chain parentChain;
    private Map<String, Object> featureMap;

    /**
     * The constructor of residues.
     * @param aminoAcid the amino acid this residue represents as Object
     * @param residueNumber a {@link Chain}-wide unique number to identify this residue
     * @see AminoAcid
     */
    public Residue(AminoAcid aminoAcid, int residueNumber) {
        this.aminoAcid = aminoAcid;
        this.residueNumber = residueNumber;
        this.atoms = new ArrayList<>();
        this.featureMap = new HashMap<>();
    }

    public Map<String, Object> getFeatureMap() {
        return featureMap;
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
     * Registers a child. This object will assign a reference to itself to the atom.
     * @param atom the atom to process
     */
    public void addAtom(Atom atom) {
        atoms.add(atom);
        atom.setParentResidue(this);
    }

    /**
     * Returns the {@link Chain}-wide unique number to identify this residue.
     * @return an integer
     */
    public int getResidueNumber() {
        return residueNumber;
    }

    /**
     * The short name of this residue. Is composed of the 1-letter code of the amino acid and the residue number. An
     * alanine at position 120 will yield <tt>A-120</tt>.
     * @return the short identifier of this residue
     * @see AminoAcid#getOneLetterCode()
     */
    public String getName() {
        return aminoAcid.getOneLetterCode() + "-" + residueNumber;
    }

    /**
     * Returns which of the 20 canonical amino acids this residue represents.
     * @return the type of this residue
     * @see AminoAcid
     */
    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }

    /**
     * Returns the {@link Chain} this residue is associated to.
     * @return the parent container
     */
    public Chain getParentChain() {
        return parentChain;
    }

    /**
     * Package-private method to set the parent reference.
     * @param parentChain the parent
     */
    void setParentChain(Chain parentChain) {
        this.parentChain = parentChain;
    }

    @Override
    public String composePDBRecord() {
        return atoms().map(Atom::composePDBRecord).collect(Collectors.joining(System.lineSeparator()));
    }

    @Override
    public String toString() {
        return this.getClass().getSimpleName() + " name='" + this.aminoAcid + "' resNum='" + this.residueNumber + "' size='" + this.atoms.size() + "'";
    }

    @Override
    public Stream<Atom> atoms() {
        return atoms.stream();
    }
}
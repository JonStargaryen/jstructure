package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

import java.util.NoSuchElementException;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The abstract representation of an amino acid.
 * max-ASA values from: http://dx.doi.org/10.1371/journal.pone.0080635
 * Created by bittrich on 5/24/17.
 */
public abstract class AminoAcid extends Group implements StandardAminoAcidIndicator {
    public static final String ALPHA_CARBON_NAME = "CA";
    public static final String BACKBONE_CARBON_NAME = "C";
    public static final String BACKBONE_NITROGEN_NAME = "N";
    public static final String BACKBONE_OXYGEN_NAME = "O";
    public static final String BACKBONE_HYDROGEN_NAME = "H";
    public static final String BETA_CARBON_NAME = "CB";
    public static final Set<String> BACKBONE_ATOM_NAMES = Stream.of(ALPHA_CARBON_NAME,
            BACKBONE_CARBON_NAME,
            BACKBONE_NITROGEN_NAME,
            BACKBONE_OXYGEN_NAME).collect(Collectors.toSet());
    /**
     * Handle to the list of names representing side getChain atoms of (any of) the standard amino acid. The intersection
     * of this 'set' and that of backbone atoms is empty.
     */
    public static final Set<String> SIDE_CHAIN_ATOM_NAMES = Stream.of("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OD1",
            "ND2", "OD2", "SG", "OE1", "NE2", "OE2", "ND1", "CD2", "CE1", "CG1", "CG2", "CD1", "CE", "NZ", "SD", "CE2",
            "OG", "OG1", "NE1", "CE3", "CZ2", "CZ3", "CH2", "OH").collect(Collectors.toSet());
    public static final Set<String> ALL_ATOM_NAMES = Stream.concat(SIDE_CHAIN_ATOM_NAMES.stream(),
            BACKBONE_ATOM_NAMES.stream()).collect(Collectors.toSet());

    /**
     * Convenience function to determine whether a certain atom is a backbone atom of an amino acid.
     * @param atom the instance to check
     * @return <code>true</code> iff this atom's name refers to that of a backbone atom
     */
    public static boolean isBackboneAtom(Atom atom) {
        String atomName = atom.getName();
        return ALPHA_CARBON_NAME.equals(atomName) ||
                BACKBONE_CARBON_NAME.equals(atomName) ||
                BACKBONE_NITROGEN_NAME.equals(atomName) ||
                BACKBONE_OXYGEN_NAME.equals(atomName) ||
                BACKBONE_HYDROGEN_NAME.equals(atomName);
    }

    public enum Family {
        ALANINE(Alanine.class, Alanine.GROUP_PROTOTYPE),
        ARGININE(Arginine.class, Arginine.GROUP_PROTOTYPE),
        ASPARAGINE(Asparagine.class, Asparagine.GROUP_PROTOTYPE),
        ASPARTIC_ACID(AsparticAcid.class, AsparticAcid.GROUP_PROTOTYPE),
        CYSTEINE(Cysteine.class, Cysteine.GROUP_PROTOTYPE),
        GLUTAMIC_ACID(GlutamicAcid.class, GlutamicAcid.GROUP_PROTOTYPE),
        GLUTAMINE(Glutamine.class, Glutamine.GROUP_PROTOTYPE),
        GLYCINE(Glycine.class, Glycine.GROUP_PROTOTYPE),
        HISTIDINE(Histidine.class, Histidine.GROUP_PROTOTYPE),
        ISOLEUCINE(Isoleucine.class, Isoleucine.GROUP_PROTOTYPE),
        LEUCINE(Leucine.class, Leucine.GROUP_PROTOTYPE),
        LYSINE(Lysine.class, Lysine.GROUP_PROTOTYPE),
        METHIONINE(Methionine.class, Methionine.GROUP_PROTOTYPE),
        PHENYLALANINE(Phenylalanine.class, Phenylalanine.GROUP_PROTOTYPE),
        PROLINE(Proline.class, Proline.GROUP_PROTOTYPE),
        SERINE(Serine.class, Serine.GROUP_PROTOTYPE),
        THREONINE(Threonine.class, Threonine.GROUP_PROTOTYPE),
        TRYPTOPHAN(Tryptophan.class, Tryptophan.GROUP_PROTOTYPE),
        TYROSINE(Tyrosine.class, Tyrosine.GROUP_PROTOTYPE),
        VALINE(Valine.class, Valine.GROUP_PROTOTYPE);

        private Class<? extends AminoAcid> representingClass;
        private GroupPrototype groupPrototype;

        Family(Class<? extends AminoAcid> representingClass, GroupPrototype groupPrototype) {
            this.representingClass = representingClass;
            this.groupPrototype = groupPrototype;
        }

        public Class<? extends AminoAcid> getRepresentingClass() {
            return representingClass;
        }

        public GroupPrototype getGroupPrototype() {
            return groupPrototype;
        }

        public static Family resolveOneLetterCode(String oneLetterCode) {
            return Stream.of(Family.values())
                    .filter(aminoAcid -> oneLetterCode.equalsIgnoreCase(aminoAcid.getGroupPrototype().getOneLetterCode().get()))
                    .findFirst()
                    .orElseThrow(() -> new NoSuchElementException("'" + oneLetterCode + "' is no valid amino acid one-letter-code"));
        }
    }

    private Atom n;
    private Atom ca;
    private Atom c;
    private Atom o;
    private Atom h;

    AminoAcid(AminoAcid aminoAcid) {
        super(aminoAcid);
        this.n = aminoAcid.n;
        this.ca = aminoAcid.ca;
        this.c = aminoAcid.c;
        this.o = aminoAcid.o;
        this.h = aminoAcid.h;
    }

    AminoAcid(String threeLetterCode,
              ResidueNumber residueNumber,
              boolean ligand) {
        super(threeLetterCode,
                residueNumber,
                ligand);
    }

    AminoAcid(String threeLetterCode,
              ResidueNumber residueNumber) {
        this(threeLetterCode, residueNumber, false);
    }

    AminoAcid(GroupPrototype groupPrototype,
              ResidueNumber residueNumber,
              boolean ligand) {
        super(groupPrototype,
                residueNumber,
                ligand);
    }

    AminoAcid(GroupPrototype groupPrototype,
              ResidueNumber residueNumber) {
        this(groupPrototype, residueNumber, false);
    }

    public Atom getN() {
        return n;
    }

    public Atom getCa() {
        return ca;
    }

    public Atom getC() {
        return c;
    }

    public Atom getO() {
        return o;
    }

    public Atom getH() {
        return h;
    }

    public String getOneLetterCode() {
        return getGroupPrototype().getOneLetterCode().orElse("?");
    }

    @Override
    protected void addAtomInternal(Atom atom) {
        if(atom.getName().equals("N") && n == null) {
            n = atom;
        }
        if(atom.getName().equals("CA") && ca == null) {
            ca = atom;
        }
        if(atom.getName().equals("C") && c == null) {
            c = atom;
        }
        if(atom.getName().equals("O") && o == null) {
            o = atom;
        }
        if(atom.getName().equals("H") && h == null) {
            h = atom;
        }
        addSideChainAtom(atom);
    }

    protected abstract void addSideChainAtom(Atom atom);

    public abstract double getMaximumAccessibleSurfaceArea();
}

package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.*;

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

    public static boolean isSideChainAtom(Atom atom) {
        return !isBackboneAtom(atom);
    }

    public enum Family {
        ALANINE(Alanine.class,
                Alanine.GROUP_PROTOTYPE,
                121.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        ARGININE(Arginine.class,
                Arginine.GROUP_PROTOTYPE,
                265.0,
                GroupPrototype.GutteridgeGrouping.GUANIDINIUM),
        ASPARAGINE(Asparagine.class,
                Asparagine.GROUP_PROTOTYPE,
                187.0,
                GroupPrototype.GutteridgeGrouping.AMIDE),
        ASPARTIC_ACID(AsparticAcid.class,
                AsparticAcid.GROUP_PROTOTYPE,
                187.0,
                GroupPrototype.GutteridgeGrouping.CARBOXYLATE),
        CYSTEINE(Cysteine.class,
                Cysteine.GROUP_PROTOTYPE,
                148.0,
                GroupPrototype.GutteridgeGrouping.THIOL),
        GLUTAMIC_ACID(GlutamicAcid.class,
                GlutamicAcid.GROUP_PROTOTYPE,
                214.0,
                GroupPrototype.GutteridgeGrouping.AMIDE),
        GLUTAMINE(Glutamine.class,
                Glutamine.GROUP_PROTOTYPE,
                214.0,
                GroupPrototype.GutteridgeGrouping.CARBOXYLATE),
        GLYCINE(Glycine.class,
                Glycine.GROUP_PROTOTYPE,
                97.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        HISTIDINE(Histidine.class,
                Histidine.GROUP_PROTOTYPE,
                216.0,
                GroupPrototype.GutteridgeGrouping.IMIDAZOLE),
        ISOLEUCINE(Isoleucine.class,
                Isoleucine.GROUP_PROTOTYPE,
                195.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        LEUCINE(Leucine.class,
                Leucine.GROUP_PROTOTYPE,
                191.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        LYSINE(Lysine.class,
                Lysine.GROUP_PROTOTYPE,
                230.0,
                GroupPrototype.GutteridgeGrouping.AMINO),
        METHIONINE(Methionine.class,
                Methionine.GROUP_PROTOTYPE,
                203.0,
                GroupPrototype.GutteridgeGrouping.THIOL),
        PHENYLALANINE(Phenylalanine.class,
                Phenylalanine.GROUP_PROTOTYPE,
                228.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        PROLINE(Proline.class,
                Proline.GROUP_PROTOTYPE,
                154.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        SERINE(Serine.class,
                Serine.GROUP_PROTOTYPE,
                143.0,
                GroupPrototype.GutteridgeGrouping.HYDROXYL),
        THREONINE(Threonine.class,
                Threonine.GROUP_PROTOTYPE,
                163.0,
                GroupPrototype.GutteridgeGrouping.HYDROXYL),
        TRYPTOPHAN(Tryptophan.class,
                Tryptophan.GROUP_PROTOTYPE,
                264.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        TYROSINE(Tyrosine.class,
                Tyrosine.GROUP_PROTOTYPE,
                255.0,
                GroupPrototype.GutteridgeGrouping.HYDROXYL),
        VALINE(Valine.class,
                Valine.GROUP_PROTOTYPE,
                165.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        UNKNOWN_AMINO_ACID(UnknownAminoAcid.class,
                UnknownAminoAcid.GROUP_PROTOTYPE,
                121.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        PYRROLYSINE(Pyrrolysine.class,
                Pyrrolysine.GROUP_PROTOTYPE,
                //TODO maximum asa value
                154.0 + 230.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        SELENOCYSTEINE(Selenocysteine.class,
                Selenocysteine.GROUP_PROTOTYPE,
                148.0,
                GroupPrototype.GutteridgeGrouping.NONE),
        SELENOMETHIONINE(Selenomethionine.class,
                Selenomethionine.GROUP_PROTOTYPE,
                203.0,
                GroupPrototype.GutteridgeGrouping.NONE);

        private Class<? extends AminoAcid> representingClass;
        private GroupPrototype groupPrototype;
        private double maximumAccessibleSurfaceArea;
        private GroupPrototype.GutteridgeGrouping gutteridgeGrouping;

        Family(Class<? extends AminoAcid> representingClass,
               GroupPrototype groupPrototype,
               double maximumAccessibleSurfaceArea,
               GroupPrototype.GutteridgeGrouping gutteridgeGrouping) {
            this.representingClass = representingClass;
            this.groupPrototype = groupPrototype;
            this.maximumAccessibleSurfaceArea = maximumAccessibleSurfaceArea;
            this.gutteridgeGrouping = gutteridgeGrouping;
        }

        public Class<? extends AminoAcid> getRepresentingClass() {
            return representingClass;
        }

        public GroupPrototype getGroupPrototype() {
            return groupPrototype;
        }

        public double getMaximumAccessibleSurfaceArea() {
            return maximumAccessibleSurfaceArea;
        }

        public GroupPrototype.GutteridgeGrouping getGutteridgeGrouping() {
            return gutteridgeGrouping;
        }

        public static Family resolveOneLetterCode(String oneLetterCode) {
            return Stream.of(Family.values())
                    .filter(aminoAcid -> oneLetterCode.equalsIgnoreCase(aminoAcid.getGroupPrototype().getOneLetterCode().get()))
                    .findFirst()
                    .orElseThrow(() -> new NoSuchElementException("'" + oneLetterCode + "' is no valid amino acid one-letter-code"));
        }

        public static Family resolveGroupPrototype(GroupPrototype groupPrototype) {
            GroupPrototype.PolymerType polymerType = groupPrototype.getPolymerType();
            if(groupPrototype.getPolymerType() != GroupPrototype.PolymerType.PEPTIDE_LINKING && polymerType !=
                    GroupPrototype.PolymerType.PEPTIDE_LIKE && polymerType !=
                    GroupPrototype.PolymerType.PEPTIDE_TERMINUS) {
                throw new UnsupportedOperationException("method only supported for amino acids - group '" +
                        groupPrototype.getThreeLetterCode() + "' has polymer type: " + polymerType);
            }

            return Stream.of(Family.values())
                    .filter(aminoAcid -> groupPrototype.equals(aminoAcid.getGroupPrototype()))
                    .findFirst()
                    .orElse(UNKNOWN_AMINO_ACID);
        }

        public static Stream<AminoAcid.Family> canonicalAminoAcids() {
            return Stream.of(values())
                    // canonical amino acids are the first 20 of the enum - TODO not the nicest way, UnknownAminoAcid is a StandardAminoAcid however
                    .limit(20);
        }
    }

    private Atom n;
    private Atom ca;
    private Atom c;
    private Atom o;
    private Atom h;

    AminoAcid(AminoAcid aminoAcid) {
        super(aminoAcid);
//        this.n = new Atom(aminoAcid.n);
//        this.ca = new Atom(aminoAcid.ca);
//        this.c = new Atom(aminoAcid.c);
//        this.o = new Atom(aminoAcid.o);
//        this.h = new Atom(aminoAcid.h);
        atoms().forEach(this::addAtomInternal);
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
}

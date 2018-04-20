package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.*;

import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The abstract representation of an amino acid.
 * max-ASA values from: http://dx.doi.org/10.1371/journal.pone.0080635
 * frequencies from: Paula Yurkanis Bruice: Organic Chemistry. Pearson Education Inc., 2004, 4. Auflage, S. 960–962, ISBN 0-13-121730-5.
 * pI: http://www.mhhe.com/physsci/chemistry/carey5e/Ch27/ch27-1-4-2.html
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

    public enum Family implements GroupFamily {
        ALANINE(Alanine.class,
                Alanine.GROUP_PROTOTYPE,
                121.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                9.0,
                6.00),
        ARGININE(Arginine.class,
                Arginine.GROUP_PROTOTYPE,
                265.0,
                GroupPrototype.GutteridgeGrouping.GUANIDINIUM,
                4.7,
                10.76),
        ASPARAGINE(Asparagine.class,
                Asparagine.GROUP_PROTOTYPE,
                187.0,
                GroupPrototype.GutteridgeGrouping.AMIDE,
                4.4,
                5.41),
        ASPARTIC_ACID(AsparticAcid.class,
                AsparticAcid.GROUP_PROTOTYPE,
                187.0,
                GroupPrototype.GutteridgeGrouping.CARBOXYLATE,
                5.5,
                2.77),
        CYSTEINE(Cysteine.class,
                Cysteine.GROUP_PROTOTYPE,
                148.0,
                GroupPrototype.GutteridgeGrouping.THIOL,
                2.8,
                5.07),
        GLUTAMIC_ACID(GlutamicAcid.class,
                GlutamicAcid.GROUP_PROTOTYPE,
                214.0,
                GroupPrototype.GutteridgeGrouping.AMIDE,
                3.9,
                3.22),
        GLUTAMINE(Glutamine.class,
                Glutamine.GROUP_PROTOTYPE,
                214.0,
                GroupPrototype.GutteridgeGrouping.CARBOXYLATE,
                6.2,
                5.65),
        GLYCINE(Glycine.class,
                Glycine.GROUP_PROTOTYPE,
                97.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                7.5,
                5.97),
        HISTIDINE(Histidine.class,
                Histidine.GROUP_PROTOTYPE,
                216.0,
                GroupPrototype.GutteridgeGrouping.IMIDAZOLE,
                2.1,
                7.59),
        ISOLEUCINE(Isoleucine.class,
                Isoleucine.GROUP_PROTOTYPE,
                195.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                4.6,
                6.02),
        LEUCINE(Leucine.class,
                Leucine.GROUP_PROTOTYPE,
                191.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                7.5,
                5.98),
        LYSINE(Lysine.class,
                Lysine.GROUP_PROTOTYPE,
                230.0,
                GroupPrototype.GutteridgeGrouping.AMINO,
                7.0,
                9.74),
        METHIONINE(Methionine.class,
                Methionine.GROUP_PROTOTYPE,
                203.0,
                GroupPrototype.GutteridgeGrouping.THIOL,
                1.7,
                5.74),
        PHENYLALANINE(Phenylalanine.class,
                Phenylalanine.GROUP_PROTOTYPE,
                228.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                3.5,
                5.48),
        PROLINE(Proline.class,
                Proline.GROUP_PROTOTYPE,
                154.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                4.6,
                6.30),
        SERINE(Serine.class,
                Serine.GROUP_PROTOTYPE,
                143.0,
                GroupPrototype.GutteridgeGrouping.HYDROXYL,
                7.1,
                5.68),
        THREONINE(Threonine.class,
                Threonine.GROUP_PROTOTYPE,
                163.0,
                GroupPrototype.GutteridgeGrouping.HYDROXYL,
                6.0,
                5.60),
        TRYPTOPHAN(Tryptophan.class,
                Tryptophan.GROUP_PROTOTYPE,
                264.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                1.1,
                5.89),
        TYROSINE(Tyrosine.class,
                Tyrosine.GROUP_PROTOTYPE,
                255.0,
                GroupPrototype.GutteridgeGrouping.HYDROXYL,
                3.5,
                5.66),
        VALINE(Valine.class,
                Valine.GROUP_PROTOTYPE,
                165.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                6.9,
                5.96),
        UNKNOWN_AMINO_ACID(UnknownAminoAcid.class,
                UnknownAminoAcid.GROUP_PROTOTYPE,
                121.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                0.0,
                0.0),
        PYRROLYSINE(Pyrrolysine.class,
                Pyrrolysine.GROUP_PROTOTYPE,
                //TODO maximum asa value
                154.0 + 230.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                0.0,
                0.0),
        SELENOCYSTEINE(Selenocysteine.class,
                Selenocysteine.GROUP_PROTOTYPE,
                148.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                0.0,
                0.0),
        SELENOMETHIONINE(Selenomethionine.class,
                Selenomethionine.GROUP_PROTOTYPE,
                203.0,
                GroupPrototype.GutteridgeGrouping.NONE,
                0.0,
                0.0);

        private Class<? extends AminoAcid> representingClass;
        private GroupPrototype groupPrototype;
        private double maximumAccessibleSurfaceArea;
        private GroupPrototype.GutteridgeGrouping gutteridgeGrouping;
        private double frequency;
        private double isoelectricPoint;

        Family(Class<? extends AminoAcid> representingClass,
               GroupPrototype groupPrototype,
               double maximumAccessibleSurfaceArea,
               GroupPrototype.GutteridgeGrouping gutteridgeGrouping,
               double frequency,
               double isoelectricPoint) {
            this.representingClass = representingClass;
            this.groupPrototype = groupPrototype;
            this.maximumAccessibleSurfaceArea = maximumAccessibleSurfaceArea;
            this.gutteridgeGrouping = gutteridgeGrouping;
            this.frequency = frequency;
            this.isoelectricPoint = isoelectricPoint;
        }

        @Override
        public Class<? extends AminoAcid> getRepresentingClass() {
            return representingClass;
        }

        @Override
        public GroupPrototype getGroupPrototype() {
            return groupPrototype;
        }

        public double getMaximumAccessibleSurfaceArea() {
            return maximumAccessibleSurfaceArea;
        }

        public GroupPrototype.GutteridgeGrouping getGutteridgeGrouping() {
            return gutteridgeGrouping;
        }

        public double getFrequency() {
            return frequency;
        }

        public double getIsoelectricPoint() {
            return isoelectricPoint;
        }

        public static AminoAcid createAminoAcid(String pdbName, ResidueIdentifier residueIdentifier, boolean ligand) {
            Class<? extends AminoAcid> representingClass = resolveThreeLetterCode(pdbName).representingClass;
            // use special constructor for UnknownAminoAcid
            if(representingClass.isAssignableFrom(UnknownAminoAcid.class)) {
                return new UnknownAminoAcid(pdbName, residueIdentifier, ligand);
            } else {
                try {
                    return representingClass.getConstructor(ResidueIdentifier.class, boolean.class).newInstance(residueIdentifier, ligand);
                } catch (Exception e) {
                    throw new RuntimeException("creation of AminoAcid instance failed", e);
                }
            }
        }

        public static Family resolveOneLetterCode(char oneLetterCode) {
            return resolveOneLetterCode(String.valueOf(oneLetterCode));
        }

        public static Family resolveOneLetterCode(String oneLetterCode) {
            return Stream.of(Family.values())
                    .filter(aminoAcid -> oneLetterCode.equalsIgnoreCase(aminoAcid.getOneLetterCode()))
                    .findFirst()
                    .orElse(Family.UNKNOWN_AMINO_ACID);
        }

        public static Family resolveThreeLetterCode(String threeLetterCode) {
            return Stream.of(Family.values())
                    .filter(aminoAcid -> threeLetterCode.equalsIgnoreCase(aminoAcid.getThreeLetterCode()))
                    .findFirst()
                    .orElse(Family.UNKNOWN_AMINO_ACID);
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
                    // canonical amino acids are the first 20 of the enum
                    .limit(20);
        }

        /**
         * Convenience method to retrieve the one-letter-code (as normally you cannot be sure of a group to have an olc,
         * but for entries of this enum we are certain.
         * @return this amino acid's one-letter-code
         */
        public String getOneLetterCode() {
            return groupPrototype.getOneLetterCode().get();
        }

        @Override
        public String getThreeLetterCode() {
            return groupPrototype.getThreeLetterCode();
        }
    }

    private Atom n;
    private Atom ca;
    private Atom c;
    private Atom o;
    private Atom h;

    public int getAminoAcidIndex() {
        // determine residue index and subtract number of non-amino acids
        int residueIndex = getResidueIndex();
        return residueIndex - (int) getParentChain().aminoAcids()
                .limit(residueIndex)
                .filter(group -> !group.isAminoAcid())
                .count();
    }

    AminoAcid(AminoAcid aminoAcid, boolean deep) {
        super(aminoAcid, deep);
        atoms().forEach(this::addAtomInternal);
    }

    AminoAcid(String threeLetterCode,
              ResidueIdentifier residueIdentifier,
              boolean ligand) {
        super(threeLetterCode,
                residueIdentifier,
                ligand);
    }

    AminoAcid(String threeLetterCode,
              ResidueIdentifier residueIdentifier) {
        this(threeLetterCode, residueIdentifier, false);
    }

    AminoAcid(GroupPrototype groupPrototype,
              ResidueIdentifier residueIdentifier,
              boolean ligand) {
        super(groupPrototype,
                residueIdentifier,
                ligand);
    }

    AminoAcid(GroupPrototype groupPrototype,
              ResidueIdentifier residueIdentifier) {
        this(groupPrototype, residueIdentifier, false);
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

    /**
     * Access to the previous amino acid in the chain.
     * <code>getPreviousAminoAcid() is equivalent to
     * </code><blockquote><pre>
     * getAminoAcidWithOffset(-1)
     * </pre></blockquote>
     * @see #getAminoAcidWithOffset(int)
     * @return an optional wrapping the previous amino acid
     */
    public Optional<AminoAcid> getPreviousAminoAcid() {
        return getAminoAcidWithOffset(-1);
    }

    /**
     * Access to the amino acid with a custom offset in the chain. Zero will return the current instance, negative
     * values will navigate towards the N-terminus, positive values to the C-terminus. Wrapped as optional as the
     * operation is by no means guaranteed to succeed, the group could not exist (as the current instance is a terminal
     * amino acid) or be no {@link AminoAcid}.
     * @see #getPreviousAminoAcid()
     * @see #getNextAminoAcid()
     * @return an optional wrapping the amino acid with the given offset
     */
    public Optional<AminoAcid> getAminoAcidWithOffset(int offset) {
        try {
            Chain chain = getParentChain();
            int index = chain.getGroups().indexOf(this);
            return Optional.of((AminoAcid) chain.getGroups().get(index + offset));
        } catch (Exception e) {
            return Optional.empty();
        }
    }

    /**
     * Access to the next amino acid in the chain.
     * <code>getPreviousAminoAcid() is equivalent to
     * </code><blockquote><pre>
     * getAminoAcidWithOffset(1)
     * </pre></blockquote>
     * @see #getAminoAcidWithOffset(int)
     * @return an optional wrapping the next amino acid
     */
    public Optional<AminoAcid> getNextAminoAcid() {
        return getAminoAcidWithOffset(1);
    }
}

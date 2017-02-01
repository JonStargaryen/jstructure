package de.bioforscher.jstructure.model.structure.family;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumerates all 20 canonical amino acids and keeps track of the present atom names for each amino acid.
 * max-ASA values from: http://dx.doi.org/10.1371/journal.pone.0080635
 * Created by S on 05.01.2017.
 */
public enum AminoAcidFamily implements AtomicFamily {
    ALANINE("ALA", "A", new String[] { "CB" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 121.0),
    ARGININE("ARG", "R", new String[] { "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2" }, AminoAcidFamily.GutteridgeGrouping.GUANIDINIUM, AminoAcidFamily.ANVILGrouping.POLAR, 265.0),
    ASPARAGINE("ASN", "N", new String[] { "CB", "CG", "OD1", "ND2" }, AminoAcidFamily.GutteridgeGrouping.AMIDE, AminoAcidFamily.ANVILGrouping.POLAR, 187.0),
    ASPARTIC_ACID("ASP", "D", new String[] { "CB", "CG", "OD1", "OD2" }, AminoAcidFamily.GutteridgeGrouping.CARBOXYLATE, AminoAcidFamily.ANVILGrouping.POLAR, 187.0),
    CYSTEINE("CYS", "C", new String[] { "CB", "SG" }, AminoAcidFamily.GutteridgeGrouping.THIOL, AminoAcidFamily.ANVILGrouping.MEMBRANE, 148.0),
    GLUTAMINE("GLN", "Q", new String[] { "CB", "CG", "CD", "OE1", "NE2" }, AminoAcidFamily.GutteridgeGrouping.AMIDE, AminoAcidFamily.ANVILGrouping.POLAR, 214.0),
    GLUTAMIC_ACID("GLU", "E", new String[] { "CB", "CG", "CD", "OE1", "OE2" }, AminoAcidFamily.GutteridgeGrouping.CARBOXYLATE, AminoAcidFamily.ANVILGrouping.POLAR, 214.0),
    GLYCINE("GLY", "G", new String[] {}, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 97.0),
    HISTIDINE("HIS", "H", new String[] { "CB", "CG", "ND1", "CD2", "CE1", "NE2" }, AminoAcidFamily.GutteridgeGrouping.IMIDAZOLE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 216.0),
    ISOLEUCINE("ILE", "I", new String[] { "CB", "CG1", "CG2", "CD1" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 195.0),
    LEUCINE("LEU", "L", new String[] { "CB", "CG", "CD1", "CD2" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 191.0),
    LYSINE("LYS", "K", new String[] { "CB", "CG", "CD", "CE", "NZ" }, AminoAcidFamily.GutteridgeGrouping.AMINO, AminoAcidFamily.ANVILGrouping.POLAR, 230.0),
    METHIONINE("MET", "M", new String[] { "CB", "CG", "SD", "CE" }, AminoAcidFamily.GutteridgeGrouping.THIOL, AminoAcidFamily.ANVILGrouping.MEMBRANE, 203.0),
    PHENYLALANINE("PHE", "F", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 228.0),
    PROLINE("PRO", "P", new String[] { "CB", "CG", "CD" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.POLAR, 154.0),
    SERINE("SER", "S", new String[] { "CB", "OG" }, AminoAcidFamily.GutteridgeGrouping.HYDROXYL, AminoAcidFamily.ANVILGrouping.MEMBRANE, 143.0),
    THREONINE("THR", "T", new String[] { "CB", "OG1", "CG2" }, AminoAcidFamily.GutteridgeGrouping.HYDROXYL, AminoAcidFamily.ANVILGrouping.MEMBRANE, 163.0),
    TRYPTOPHAN("TRP", "W", new String[] { "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.POLAR, 264.0),
    TYROSINE("TYR", "Y", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH" }, AminoAcidFamily.GutteridgeGrouping.HYDROXYL, AminoAcidFamily.ANVILGrouping.POLAR, 255.0),
    VALINE("VAL", "V", new String[] { "CB", "CG1", "CG2" }, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.MEMBRANE, 165.0),
    UNKNOWN("UNK", "X", new String[] {}, AminoAcidFamily.GutteridgeGrouping.NONE, AminoAcidFamily.ANVILGrouping.UNKNOWN, 121.0);

    static Logger logger = LoggerFactory.getLogger(AminoAcidFamily.class);

    public boolean isStandardAminoAcid() {
        return !oneLetterCode.equals(UNKNOWN.getOneLetterCode());
    }

    public double getMaximalAccessibleSurfaceArea() {
        return maximalAccessibleSurfaceArea;
    }

    /**
     * Gutteridge grouping (doi:10.1016/j.tibs.2005.09.006)
     */
    public enum GutteridgeGrouping {
        NONE,
        THIOL,
        HYDROXYL,
        GUANIDINIUM,
        IMIDAZOLE,
        AMINO,
        CARBOXYLATE,
        AMIDE
    }

    /**
     * The grouping by the ANVIL algorithm according to the tendency of amino acids to be placed within the membrane layer.
     */
    public enum ANVILGrouping {
        MEMBRANE,
        POLAR,
        UNKNOWN
    }

    public interface ATOM_NAMES {
        /**
         * Handle to the list of names representing hydrogen atoms.
         */
        Set<String> H_ATOM_NAMES = Stream.of("H", "D", "T").collect(Collectors.toSet());

        /**
         * Handle to the list of names representing alpha carbons.
         */
        String CA_ATOM_NAME = "CA";
        String C_ATOM_NAME = "C";
        String N_ATOM_NAME = "N";
        String O_ATOM_NAME = "O";

        /**
         * Handle to the list of names representing backbone atoms.
         */
        Set<String> BACKBONE_ATOM_NAMES = Stream.of("N", "CA", "C", "O").collect(Collectors.toSet());

        /**
         * The beta carbon name.
         */
        String CB_ATOM_NAME = "CB";

        /**
         * Handle to the list of names representing side getChain atoms of (any of) the standard amino acid. The intersection
         * of this 'set' and that of backbone atoms is empty.
         */
        Set<String> SIDECHAIN_ATOM_NAMES = Stream.of("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OD1", "ND2",
                "OD2", "SG", "OE1", "NE2", "OE2", "ND1", "CD2", "CE1", "CG1", "CG2", "CD1", "CE", "NZ", "SD", "CE2",
                "OG", "OG1", "NE1", "CE3", "CZ2", "CZ3", "CH2", "OH").collect(Collectors.toSet());

        /**
         * All occurring atom names (ignoring hydrogen atoms).
         */
        Set<String> ALL_ATOM_NAMES = Stream.concat(SIDECHAIN_ATOM_NAMES.stream(),
                BACKBONE_ATOM_NAMES.stream()).collect(Collectors.toSet());
    }

    private static Map<String, AminoAcidFamily> allAminoAcids;

    static {
        /*
         * register all amino acids with their one-/three-letter code as well as full name, so the proper amino acid can
         * be retrieved easily by a String
         */
        allAminoAcids = new HashMap<>();
        for (AminoAcidFamily aa : AminoAcidFamily.values()){
            allAminoAcids.put(aa.oneLetterCode.toLowerCase(), aa);
            allAminoAcids.put(aa.threeLetterCode.toLowerCase(), aa);
            allAminoAcids.put(aa.toString().toLowerCase(), aa);
        }
    }

    private final String oneLetterCode;
    private final String threeLetterCode;
    private final List<String> sideChainAtomNames;
    private final AminoAcidFamily.GutteridgeGrouping gutteridgeGrouping;
    private final AminoAcidFamily.ANVILGrouping anvilGrouping;
    private final AminoAcidFamily.AtomNameComparator atomNameComparator;
    private final double maximalAccessibleSurfaceArea;

    AminoAcidFamily(String threeLetterCode, String oneLetterCode, String[] sideChainAtomNames,
              AminoAcidFamily.GutteridgeGrouping gutteridgeGrouping, AminoAcidFamily.ANVILGrouping anvilGrouping, double maximalAccessibleSurfaceArea) {
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
        this.sideChainAtomNames = Arrays.asList(sideChainAtomNames);
        this.atomNameComparator = new AminoAcidFamily.AtomNameComparator(allAtomNames().collect(Collectors.toList()));
        this.gutteridgeGrouping = gutteridgeGrouping;
        this.anvilGrouping = anvilGrouping;
        this.maximalAccessibleSurfaceArea = maximalAccessibleSurfaceArea;
    }

    /**
     * Convenience access to {@link AminoAcidFamily} objects based on their name.
     * @param aminoAcidName the one-/three-letter code or full name of a amino acid
     * @return the proper amino acid instance
     * @throws IllegalArgumentException if the argument does not correspond to any of the 20 canonical amino acids
     */
    @Deprecated
    public static Optional<AminoAcidFamily> valueOfIgnoreCase(String aminoAcidName) throws IllegalArgumentException {
        return Optional.ofNullable(allAminoAcids.get(aminoAcidName.toLowerCase()));
    }

    /**
     * Access to the one-letter code of this amino acid.
     * @return the proper one-letter code, e.g. "A"
     */
    public String getOneLetterCode() {
        return oneLetterCode;
    }

    /**
     * Access to the three-letter code of this amino acid.
     * @return the proper three-letter code, e.g. "ALA"
     */
    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    /**
     * Access to the side getChain atom names of this amino acid.
     * @return a select containing all side getChain atom names - that of {@link AminoAcidFamily#GLYCINE} is empty
     */
    public Stream<String> sideChainAtomNames() {
        return getSideChainAtomNames().stream();
    }

    /**
     * Access to the side getChain atom names of this amino acid.
     * @return a list containing all side getChain atom names - that of {@link AminoAcidFamily#GLYCINE} is empty
     */
    public List<String> getSideChainAtomNames() {
        return sideChainAtomNames;
    }

    /**
     * Access to all atom names present in this amino acid.
     * @return a list containing all atom names of this amino acid
     */
    public Stream<String> allAtomNames() {
        return Stream.concat(AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES.stream(), sideChainAtomNames.stream());
    }

    /**
     * Grouping of amino acids by their functional group.
     * doi:10.1016/j.tibs.2005.09.006
     * @return the enum representing this group
     */
    public AminoAcidFamily.GutteridgeGrouping getGutteridgeGrouping() {
        return gutteridgeGrouping;
    }

    /**
     * Grouping of amino acids by their tendency to be embedded in the membrane.
     * @return the enum representing this group
     */
    public AminoAcidFamily.ANVILGrouping getANVILGrouping() {
        return anvilGrouping;
    }

    /**
     * Orders all atoms associated to this {@link AtomContainer} in an reproducible way in agreement with the ordering
     * in <tt>PDB</tt> files.
     */
    public static Group orderAtomsByName(Group group) {
        if(!group.isAminoAcid()) {
            throw new IllegalArgumentException("cannot sort atoms by name for non-amino acids");
        }
        group.getAtoms().sort(AminoAcidFamily.valueOfIgnoreCase(group.getThreeLetterCode()).get().atomNameComparator);
        return group;
    }

    /**
     * Can order atoms within a container according the ordering in <tt>PDB</tt> files.
     */
    final static class AtomNameComparator implements Comparator<Atom> {
        private final String[] atomNames;

        AtomNameComparator(List<String> allAminoAcidAtomNames) {
            atomNames = allAminoAcidAtomNames.toArray(new String[allAminoAcidAtomNames.size()]);
        }

        @Override
        public int compare(Atom a1, Atom a2) {
            return orderOf(a1) - orderOf(a2);
        }

        private int orderOf(Atom atom) {
            for(int i = 0; i < atomNames.length; i++) {
                if(atomNames[i].equals(atom.getName())) {
                    return i;
                }
            }
            return atomNames.length;
        }
    }
}

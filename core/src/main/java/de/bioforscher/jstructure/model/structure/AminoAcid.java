package de.bioforscher.jstructure.model.structure;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumerates all 20 canonical amino acids and keeps track of the present atom names for each amino acid.
 * Created by S on 28.09.2016.
 */
public enum AminoAcid {
    ALANINE("ALA", "A", new String[] { "CB" }, GutteridgeGrouping.NONE, ANVILGrouping.MEMBRANE),
	ARGININE("ARG", "R", new String[] { "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2" }, GutteridgeGrouping.GUANIDINIUM, ANVILGrouping.POLAR),
	ASPARAGINE("ASN", "N", new String[] { "CB", "CG", "OD1", "ND2" }, GutteridgeGrouping.AMIDE, ANVILGrouping.POLAR),
    ASPARTIC_ACID("ASP", "D", new String[] { "CB", "CG", "OD1", "OD2" }, GutteridgeGrouping.CARBOXYLATE, ANVILGrouping.POLAR),
	CYSTEINE("CYS", "C", new String[] { "CB", "SG" }, GutteridgeGrouping.THIOL, ANVILGrouping.MEMBRANE),
	GLUTAMINE("GLN", "Q", new String[] { "CB", "CG", "CD", "OE1", "NE2" }, GutteridgeGrouping.AMIDE, ANVILGrouping.POLAR),
	GLUTAMIC_ACID("GLU", "E", new String[] { "CB", "CG", "CD", "OE1", "OE2" }, GutteridgeGrouping.CARBOXYLATE, ANVILGrouping.POLAR),
	GLYCINE("GLY", "G", new String[] {}, GutteridgeGrouping.NONE, ANVILGrouping.MEMBRANE),
	HISTIDINE("HIS", "H", new String[] { "CB", "CG", "ND1", "CD2", "CE1", "NE2" }, GutteridgeGrouping.IMIDAZOLE, ANVILGrouping.MEMBRANE),
	ISOLEUCINE("ILE", "I", new String[] { "CB", "CG1", "CG2", "CD1" }, GutteridgeGrouping.NONE, ANVILGrouping.MEMBRANE),
	LEUCINE("LEU", "L", new String[] { "CB", "CG", "CD1", "CD2" }, GutteridgeGrouping.NONE, ANVILGrouping.MEMBRANE),
	LYSINE("LYS", "K", new String[] { "CB", "CG", "CD", "CE", "NZ" }, GutteridgeGrouping.AMINO, ANVILGrouping.POLAR),
	METHIONINE("MET", "M", new String[] { "CB", "CG", "SD", "CE" }, GutteridgeGrouping.THIOL, ANVILGrouping.MEMBRANE),
	PHENYLALANINE("PHE", "F", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ" }, GutteridgeGrouping.NONE, ANVILGrouping.MEMBRANE),
	PROLINE("PRO", "P", new String[] { "CB", "CG", "CD" }, GutteridgeGrouping.NONE, ANVILGrouping.POLAR),
	SERINE("SER", "S", new String[] { "CB", "OG" }, GutteridgeGrouping.HYDROXYL, ANVILGrouping.MEMBRANE),
	THREONINE("THR", "T", new String[] { "CB", "OG1", "CG2" }, GutteridgeGrouping.HYDROXYL, ANVILGrouping.MEMBRANE),
	TRYPTOPHAN("TRP", "W", new String[] { "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2" }, GutteridgeGrouping.NONE, ANVILGrouping.POLAR),
	TYROSINE("TYR", "Y", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH" }, GutteridgeGrouping.HYDROXYL, ANVILGrouping.POLAR),
    VALINE("VAL", "V", new String[] { "CB", "CG1", "CG2" }, GutteridgeGrouping.NONE, ANVILGrouping.MEMBRANE),
    UNKNOWN("UNK", "X", new String[] {}, GutteridgeGrouping.NONE, ANVILGrouping.UNKNOWN);

    static Logger logger = LoggerFactory.getLogger(AminoAcid.class);

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
        List<String> H_ATOM_NAME = Arrays.asList("H", "D", "T");
        /**
         * Handle to the list of names representing alpha carbons.
         */
        List<String> CA_ATOM_NAMES = Collections.singletonList("CA");
        List<String> C_ATOM_NAMES = Collections.singletonList("C");
        List<String> N_ATOM_NAMES = Collections.singletonList("N");
        List<String> O_ATOM_NAMES = Collections.singletonList("O");
        /**
         * Handle to the list of names representing backbone atoms.
         */
        List<String> BACKBONE_ATOM_NAMES = Arrays.asList("N", "CA", "C", "O");
        /**
         * Handle to the list of names representing side chain atoms of (any of) the standard amino acid. The intersection
         * of this 'set' and that of backbone atoms is empty.
         */
        List<String> SIDECHAIN_ATOM_NAMES = Arrays.asList("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OD1", "ND2",
                "OD2", "SG", "OE1", "NE2", "OE2", "ND1", "CD2", "CE1", "CG1", "CG2", "CD1", "CE", "NZ", "SD", "CE2",
                "OG", "OG1", "NE1", "CE3", "CZ2", "CZ3", "CH2", "OH"
        );
    }

    private static Map<String, AminoAcid> allAminoAcids;

    static {
        /*
         * register all amino acids with their one-/three-letter code as well as full name, so the proper amino acid can
         * be retrieved easily by a String
         */
        allAminoAcids = new HashMap<>();
        for (AminoAcid aa : AminoAcid.values()){
            allAminoAcids.put(aa.oneLetterCode.toLowerCase(), aa);
            allAminoAcids.put(aa.threeLetterCode.toLowerCase(), aa);
            allAminoAcids.put(aa.toString().toLowerCase(), aa);
        }
    }

    private final String oneLetterCode;
    private final String threeLetterCode;
    private final List<String> sideChainAtomNames;
    private final GutteridgeGrouping gutteridgeGrouping;
    private final ANVILGrouping anvilGrouping;
    private final AtomNameComparator atomNameComparator;

    AminoAcid(String threeLetterCode, String oneLetterCode, String[] sideChainAtomNames,
              GutteridgeGrouping gutteridgeGrouping, ANVILGrouping anvilGrouping) {
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
        this.sideChainAtomNames = Arrays.asList(sideChainAtomNames);
        this.atomNameComparator = new AtomNameComparator(allAtomNames().collect(Collectors.toList()));
        this.gutteridgeGrouping = gutteridgeGrouping;
        this.anvilGrouping = anvilGrouping;
    }

    /**
     * Convenience access to {@link AminoAcid} objects based on their name.
     * @param aminoAcidName the one-/three-letter code or full name of a amino acid
     * @return the proper amino acid instance
     * @throws IllegalArgumentException if the argument does not correspond to any of the 20 canonical amino acids
     */
    public static AminoAcid valueOfIgnoreCase(String aminoAcidName) throws IllegalArgumentException {
        AminoAcid aa = allAminoAcids.get(aminoAcidName.toLowerCase());
        if (aa != null) {
            return aa;
        }
        //TODO handling of non-standard amino acids
        logger.warn("encountered non-standard amino acid {} - falling back to {}", aminoAcidName, AminoAcid.UNKNOWN);
//        throw new IllegalArgumentException("Invalid amino acid name: " + aminoAcidName);
        return AminoAcid.UNKNOWN;
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
     * Access to the side chain atom names of this amino acid.
     * @return a stream containing all side chain atom names - that of {@link AminoAcid#GLYCINE} is empty
     */
    public Stream<String> sideChainAtomNames() {
        return getSideChainAtomNames().stream();
    }

    /**
     * Access to the side chain atom names of this amino acid.
     * @return a list containing all side chain atom names - that of {@link AminoAcid#GLYCINE} is empty
     */
    public List<String> getSideChainAtomNames() {
        return sideChainAtomNames;
    }

    /**
     * Access to all atom names present in this amino acid.
     * @return a list containing all atom names of this amino acid
     */
    public Stream<String> allAtomNames() {
        return Stream.concat(ATOM_NAMES.BACKBONE_ATOM_NAMES.stream(), sideChainAtomNames.stream());
    }

    /**
     * Grouping of amino acids by their functional group.
     * doi:10.1016/j.tibs.2005.09.006
     * @return the enum representing this group
     */
    public GutteridgeGrouping getGutteridgeGrouping() {
        return gutteridgeGrouping;
    }

    /**
     * Grouping of amino acids by their tendency to be embedded in the membrane.
     * @return the enum representing this group
     */
    public ANVILGrouping getANVILGrouping() {
        return anvilGrouping;
    }

    /**
     * Orders all atoms associated to this {@link AtomContainer} in an reproducible way in agreement with the ordering
     * in <tt>PDB</tt> files.
     */
    public static void orderAtomsByName(Residue residue) {
        //TODO implement
//        Collections.sort(residue.atoms, aminoAcid.atomNameComparator);
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
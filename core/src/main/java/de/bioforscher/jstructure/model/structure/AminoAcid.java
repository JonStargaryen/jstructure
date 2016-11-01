package de.bioforscher.jstructure.model.structure;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumerates all 20 canonical amino acids and keeps track of the present atom names for each amino acid.
 * Created by S on 28.09.2016.
 */
public enum AminoAcid {
    ALANINE("ALA", "A", new String[] { "CB" }),
	ARGININE("ARG", "R", new String[] { "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2" }),
	ASPARAGINE("ASN", "N", new String[] { "CB", "CG", "OD1", "ND2" }),
    ASPARTIC_ACID("ASP", "D", new String[] { "CB", "CG", "OD1", "OD2" }),
	CYSTEINE("CYS", "C", new String[] { "CB", "SG" }),
	GLUTAMINE("GLN", "Q", new String[] { "CB", "CG", "CD", "OE1", "NE2" }),
	GLUTAMIC_ACID("GLU", "E", new String[] { "CB", "CG", "CD", "OE1", "OE2" }),
	GLYCINE("GLY", "G", new String[] {}),
	HISTIDINE("HIS", "H", new String[] { "CB", "CG", "ND1", "CD2", "CE1", "NE2" }),
	ISOLEUCINE("ILE", "I", new String[] { "CB", "CG1", "CG2", "CD1" }),
	LEUCINE("LEU", "L", new String[] { "CB", "CG", "CD1", "CD2" }),
	LYSINE("LYS", "K", new String[] { "CB", "CG", "CD", "CE", "NZ" }),
	METHIONINE("MET", "M", new String[] { "CB", "CG", "SD", "CE" }),
	PHENYLALANINE("PHE", "F", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ" }),
	PROLINE("PRO", "P", new String[] { "CB", "CG", "CD" }),
	SERINE("SER", "S", new String[] { "CB", "OG" }),
	THREONINE("THR", "T", new String[] { "CB", "OG1", "CG2" }),
	TRYPTOPHAN("TRP", "W", new String[] { "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2" }),
	TYROSINE("TYR", "Y", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH" }),
    VALINE("VAL", "V", new String[] { "CB", "CG1", "CG2" }),
    UNKNOWN("UNK", "X", new String[] {});

    private static Map<String, AminoAcid> allAminoAcids;

    public interface ATOM_NAMES {
        /**
         * Handle to the list of names representing hydrogen atoms.
         */
        List<String> H_ATOM_NAME = Arrays.asList("H", "D", "T");
        /**
         * Handle to the list of names representing alpha carbons.
         */
        List<String> CA_ATOM_NAMES = Arrays.asList("CA");
        List<String> C_ATOM_NAMES = Arrays.asList("C");
        List<String> N_ATOM_NAMES = Arrays.asList("N");
        List<String> O_ATOM_NAMES = Arrays.asList("O");
        /**
         * Handle to the list of names representing backbone atoms.
         */
        //TODO check whether side-chain atoms are also labeled as H
        List<String> BACKBONE_ATOM_NAMES = Arrays.asList("N", "CA", "C", "O", "H");
        /**
         * Handle to the list of names representing side chain atoms of (any of) the standard amino acid. The intersection
         * of this 'set' and that of backbone atoms is empty.
         */
        List<String> SIDECHAIN_ATOM_NAMES = Arrays.asList("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OD1", "ND2",
                "OD2", "SG", "OE1", "NE2", "OE2", "ND1", "CD2", "CE1", "CG1", "CG2", "CD1", "CE", "NZ", "SD", "CE2",
                "OG", "OG1", "NE1", "CE3", "CZ2", "CZ3", "CH2", "OH"
        );
    }

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
    private final AtomNameComparator atomNameComparator;

    AminoAcid(String threeLetterCode, String oneLetterCode, String[] sideChainAtomNames) {
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
        this.sideChainAtomNames = Arrays.asList(sideChainAtomNames);
        this.atomNameComparator = new AtomNameComparator(allAtomNames().collect(Collectors.toList()));
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
        //TODO warnings?
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
     * @return a list containing all side chain atom names - that of {@link AminoAcid#GLYCINE} is empty
     */
    public Stream<String> sideChainAtomNames() {
        return sideChainAtomNames.stream();
    }

    /**
     * Access to all atom names present in this amino acid.
     * @return a list containing all atom names of this amino acid
     */
    public Stream<String> allAtomNames() {
        return Stream.concat(ATOM_NAMES.BACKBONE_ATOM_NAMES.stream(), sideChainAtomNames.stream());
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
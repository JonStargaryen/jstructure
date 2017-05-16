package de.bioforscher.jstructure.model.structure.family;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

/**
 * The collection of nucleotides.
 * Created by S on 05.01.2017.
 */
public enum NucleotideFamily implements AtomicFamily {
    ADENOSINE("A","A"),
    DESOXYADENOSINE("A", "dA"),
    GUANOSINE("G","G"),
    DESOXYGUANOSINE("G", "dG"),
    THYMIDINE("T","T"),
    DESOXYTHYMIDINE("T", "dT"),
    URIDINE("U","U"),
    DESOXYURIDINE("U", "dU"),
    CYTIDINE("C","C"),
    DESOXYCYTIDINE("C", "dC"),
    MODIFIED_NUCLEOTIDE("M", "M"),
    UNKNOWN("X", "XX");

    private static Map<String, NucleotideFamily> allNucleotides;

    static {
        /*
         * register all amino acids with their one-/three-letter code as well as full name, so the proper amino acid can
         * be retrieved easily by a String
         */
        allNucleotides = new HashMap<>();
        for (NucleotideFamily nucleotide : NucleotideFamily.values()){
            allNucleotides.put(nucleotide.oneLetterCode.toLowerCase(), nucleotide);
            allNucleotides.put(nucleotide.threeLetterCode.toLowerCase(), nucleotide);
            allNucleotides.put(nucleotide.toString().toLowerCase(), nucleotide);
        }
    }

    private String oneLetterCode;
    private String threeLetterCode;

    NucleotideFamily(String oneLetterCode, String threeLetterCode) {
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
    }

    @Override
    public String getOneLetterCode() {
        return this.oneLetterCode;
    }

    @Override
    public String getThreeLetterCode() {
        return this.threeLetterCode;
    }

    static Optional<NucleotideFamily> valueOfIgnoreCase(String nucleotideName) {
        return Optional.ofNullable(allNucleotides.get(nucleotideName));
    }
}

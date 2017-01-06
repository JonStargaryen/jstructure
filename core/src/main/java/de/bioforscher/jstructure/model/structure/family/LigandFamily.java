package de.bioforscher.jstructure.model.structure.family;

/**
 * The rather empty collection of ligands.
 * Created by S on 05.01.2017.
 */
public enum LigandFamily implements AtomicFamily {
    UNKNOWN("*", "HET");

    private String oneLetterCode;
    private String threeLetterCode;

    LigandFamily(String oneLetterCode, String threeLetterCode) {
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
}

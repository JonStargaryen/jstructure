package de.bioforscher.jstructure.model.structure;

/**
 * Created by S on 25.05.2017.
 */
public class ResidueNumber {
    private final int residueNumber;
    private final String insertionCode;

    public ResidueNumber(int residueNumber) {
        this(residueNumber,
                "");
    }

    public ResidueNumber(int residueNumber,
                         String insertionCode) {
        this.residueNumber = residueNumber;
        this.insertionCode = insertionCode;
    }

    public int getResidueNumber() {
        return residueNumber;
    }

    public String getInsertionCode() {
        return insertionCode;
    }

    public boolean hasInsertionCode() {
        return !insertionCode.isEmpty();
    }
}

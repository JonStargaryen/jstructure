package de.bioforscher.jstructure.model.structure.identifier;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Represents the numbering of {@link Group} instances within a protein chain.
 * TODO reference to parent ChainIdentifer?
 * TODO support for residue ranges?
 * Created by S on 25.05.2017.
 */
public class ResidueIdentifier {
    private int residueNumber;
    private String insertionCode;

    ResidueIdentifier(int residueNumber) {
        this(residueNumber,
                "");
    }

    ResidueIdentifier(int residueNumber, String insertionCode) {
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

    @Override
    public String toString() {
        return residueNumber + insertionCode;
    }
}

package de.bioforscher.jstructure.model.structure.identifier;

import de.bioforscher.jstructure.model.structure.Group;

import java.util.Comparator;

/**
 * Represents the numbering of {@link Group} instances within a protein chain. Instances of this class are immutable.
 * TODO support for residue ranges?
 * Created by S on 25.05.2017.
 */
public class ResidueIdentifier {
    public static final Comparator<? super ResidueIdentifier> RESIDUE_IDENTIFIER_COMPARATOR = Comparator.comparingInt(ResidueIdentifier::getResidueNumber);
    public static final Comparator<? super Group> GROUP_COMPARATOR = Comparator.comparingInt(group -> group.getResidueIdentifier().getResidueNumber());
    private final int residueNumber;
    private final String insertionCode;

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

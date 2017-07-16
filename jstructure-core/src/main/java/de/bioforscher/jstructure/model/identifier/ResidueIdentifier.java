package de.bioforscher.jstructure.model.identifier;

import de.bioforscher.jstructure.model.structure.Group;

import java.util.Comparator;

/**
 * Represents the numbering of {@link Group} instances within a protein chain. Instances of this class are immutable.
 * Created by S on 25.05.2017.
 */
public class ResidueIdentifier extends AbstractIdentifier {
    private final int residueNumber;
    private final String insertionCode;
    public static final Comparator<? super ResidueIdentifier> RESIDUE_IDENTIFIER_COMPARATOR =
            Comparator.comparingInt(ResidueIdentifier::getResidueNumber)
                    .thenComparing(ResidueIdentifier::getInsertionCode);
    public static final Comparator<? super Group> GROUP_COMPARATOR =
            (g1, g2) -> RESIDUE_IDENTIFIER_COMPARATOR.compare(g1.getResidueIdentifier(), g2.getResidueIdentifier());

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

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;

        ResidueIdentifier that = (ResidueIdentifier) other;

        if (residueNumber != that.residueNumber) return false;
        return insertionCode.equals(that.insertionCode);
    }

    @Override
    public int hashCode() {
        int result = residueNumber;
        result = 31 * result + insertionCode.hashCode();
        return result;
    }
}

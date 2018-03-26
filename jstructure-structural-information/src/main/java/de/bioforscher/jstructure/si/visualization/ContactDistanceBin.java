package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Optional;
import java.util.stream.Stream;

/**
 * Represents the 3 states of long-range contacts in CASP experiments.
 */
public enum ContactDistanceBin {
// https://onlinelibrary.wiley.com/doi/full/10.1002/prot.24340
// short‐, medium‐ and long‐range contacts are between residues separated by 6 to 11, 12 to 23, and at least 24
    SHORT(6, 11),
    MEDIUM(12, 23),
    LONG(24, Integer.MAX_VALUE);

    private final int lowerBound;
    private final int upperBound;

    ContactDistanceBin(int lowerBound, int upperBound) {
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
    }

    public int getLowerBound() {
        return lowerBound;
    }

    public int getUpperBound() {
        return upperBound;
    }

    public boolean matchesCriterion(Pair<AminoAcid, AminoAcid> contact) {
        return matchesCriterion(contact.getLeft().getResidueIndex(), contact.getRight().getResidueIndex());
    }

    public boolean matchesCriterion(int residueIndex1, int residueIndex2) {
        int separation = Math.abs(residueIndex1 - residueIndex2);
        return separation >= lowerBound && separation <= upperBound;
    }

    public static Optional<ContactDistanceBin> resolve(Pair<ResidueIdentifier, ResidueIdentifier> contact) {
        return resolve(contact.getLeft().getResidueNumber(), contact.getRight().getResidueNumber());
    }

    public static Optional<ContactDistanceBin> resolve(int residueIndex1, int residueIndex2) {
        return Stream.of(ContactDistanceBin.values())
                .filter(contactDistanceBin -> contactDistanceBin.matchesCriterion(residueIndex1, residueIndex2))
                .findFirst();
    }
}

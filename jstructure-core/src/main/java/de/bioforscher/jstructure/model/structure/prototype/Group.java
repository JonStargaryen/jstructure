package de.bioforscher.jstructure.model.structure.prototype;

/**
 * The abstract representation of a group within a molecular structure.
 * Created by bittrich on 5/24/17.
 */
public interface Group {
    enum PolymerType {
        // non-polymer ligands
        NON_POLYMER,
        // amino acids
        PEPTIDE_LINKING,
        // nucleotides
        NA_LINKING
    }

    String getThreeLetterCode();

    Group getPrototype();

    default PolymerType getPolymerType() {
        return getPrototype().getPolymerType();
    }

    default boolean isAminoAcid() {
        return !isLigand() && getPolymerType() == PolymerType.PEPTIDE_LINKING;
    }

    default boolean isNucleotide() {
        return !isLigand() && getPolymerType() == PolymerType.NA_LINKING;
    }

    default boolean isWater() {
        return getThreeLetterCode().equals("HOH");
    }

    boolean isLigand();
}

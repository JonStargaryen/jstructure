package de.bioforscher.jstructure.model.structure.nucleotide;

/**
 * Marks non-standard nucleotides.
 * Created by bittrich on 5/30/17.
 */
public interface NonStandardNucleotide extends StandardNucleotideIndicator {
    @Override
    default boolean isStandardNucleotide() {
        return false;
    }

    Class<? extends StandardNucleotide> getParentNucleotide();
}

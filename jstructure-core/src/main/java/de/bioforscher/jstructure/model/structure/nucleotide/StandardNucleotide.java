package de.bioforscher.jstructure.model.structure.nucleotide;

/**
 * Marks standard nucleotides.
 * Created by bittrich on 5/30/17.
 */
public interface StandardNucleotide extends StandardNucleotideIndicator {
    @Override
    default boolean isStandardNucleotide() {
        return true;
    }
}

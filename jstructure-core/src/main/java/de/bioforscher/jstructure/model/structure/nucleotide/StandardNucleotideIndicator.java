package de.bioforscher.jstructure.model.structure.nucleotide;

/**
 * Reports whether this nucleotide is standard or not.
 * Created by bittrich on 5/30/17.
 */
public interface StandardNucleotideIndicator {
    /**
     * True iff this nucleotide is one of the 8 standard ones.
     * @return <code>true</code> for standard nucleotides
     */
    boolean isStandardNucleotide();
}

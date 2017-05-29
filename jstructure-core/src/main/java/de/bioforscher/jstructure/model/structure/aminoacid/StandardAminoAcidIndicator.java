package de.bioforscher.jstructure.model.structure.aminoacid;

/**
 * Reports whether this amino acid is one of the 20 canonical or not.
 * Created by bittrich on 5/24/17.
 */
public interface StandardAminoAcidIndicator {
    /**
     * True iff this amino acid is one of the 20 canonical.
     * @return <code>true</code> for canonical amino acids
     */
    boolean isStandardAminoAcid();
}

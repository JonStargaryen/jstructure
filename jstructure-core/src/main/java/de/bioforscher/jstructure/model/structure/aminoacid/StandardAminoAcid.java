package de.bioforscher.jstructure.model.structure.aminoacid;

/**
 * Marks standard amino acids.
 * Created by bittrich on 5/24/17.
 */
public interface StandardAminoAcid extends StandardAminoAcidIndicator {
    @Override
    default boolean isStandardAminoAcid() {
        return true;
    }
}
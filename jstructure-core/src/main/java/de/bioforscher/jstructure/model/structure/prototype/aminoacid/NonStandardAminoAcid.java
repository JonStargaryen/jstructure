package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

/**
 * Marks non-standard amino acids.
 * Created by bittrich on 5/24/17.
 */
public interface NonStandardAminoAcid extends StandardAminoAcidIndicator {
    @Override
    default boolean isStandardAminoAcid() {
        return false;
    }

    Class<? extends StandardAminoAcid> getParentAminoAcid();
}

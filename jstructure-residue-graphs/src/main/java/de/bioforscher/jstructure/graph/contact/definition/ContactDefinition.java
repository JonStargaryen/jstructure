package de.bioforscher.jstructure.graph.contact.definition;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public interface ContactDefinition {
    boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2);

    default String composeCaspRRLine(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        throw new UnsupportedOperationException(this.getClass().getSimpleName() + " does not support Confold reconstruction");
    }

    default String getConfoldRRType() {
        throw new UnsupportedOperationException(this.getClass().getSimpleName() + " does not support Confold reconstruction");
    }
}

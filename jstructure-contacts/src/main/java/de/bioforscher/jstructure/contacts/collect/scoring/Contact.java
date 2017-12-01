package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class Contact extends Edge<AminoAcid> {
    private final String type;

    public Contact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        this(aminoAcid1, aminoAcid2, 1.0);
    }

    public Contact(AminoAcid aminoAcid1, AminoAcid aminoAcid2, double weight) {
        this(aminoAcid1, aminoAcid2, weight, "undefined");
    }

    public Contact(AminoAcid aminoAcid1, AminoAcid aminoAcid2, String type) {
        this(aminoAcid1, aminoAcid2, 1.0, type);
    }

    public Contact(AminoAcid aminoAcid1, AminoAcid aminoAcid2, double weight, String type) {
        super(aminoAcid1, aminoAcid2, weight);
        this.type = type;
    }

    public String getType() {
        return type;
    }
}

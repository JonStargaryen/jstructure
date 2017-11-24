package de.bioforscher.jstructure.contacts.visualization;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class JsonNode {
    private final String name;

    public JsonNode(AminoAcid aminoAcid) {
        this.name = aminoAcid.getResidueIdentifier().getResidueNumber() + "-" + aminoAcid.getOneLetterCode();
    }

    public String getName() {
        return name;
    }
}

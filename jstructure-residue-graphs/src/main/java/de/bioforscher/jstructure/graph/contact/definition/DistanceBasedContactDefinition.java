package de.bioforscher.jstructure.graph.contact.definition;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public abstract class DistanceBasedContactDefinition implements ContactDefinition {
    private final double distance;
    protected final double squaredDistance;

    DistanceBasedContactDefinition(double distance) {
        this.distance = distance;
        this.squaredDistance = distance * distance;
    }

    public double getDistance() {
        return distance;
    }

    @Override
    public String composeCaspRRLine(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return (aminoAcid1.getAminoAcidIndex() + 1) + " " + (aminoAcid2.getAminoAcidIndex() + 1) + " 0 " + StandardFormat.format(distance) + " 1";
    }
}

package de.bioforscher.jstructure.feature.sse;

/**
 * The enum of secondary structure elements.
 * Created by S on 04.10.2016.
 */
public enum SecondaryStructureType {
    COIL("c"),
    BEND("S"),
    TURN("T"),
    POLYPROLINE_HELIX("P"),
    PI_HELIX("I"),
    THREE_TEN_HELIX("G"),
    BRIDGE("B"),
    EXTENDED("E"),
    ALPHA_HELIX("H");

    private String oneLetterRepresentation;

    SecondaryStructureType(String oneLetterRepresentation) {
        this.oneLetterRepresentation = oneLetterRepresentation;
    }

    public String getOneLetterRepresentation() {
        return oneLetterRepresentation;
    }

    public String getReducedRepresentation() {
        if(isHelixType()) {
            return ALPHA_HELIX.oneLetterRepresentation;
        }

        if(isStrandType()) {
            return EXTENDED.oneLetterRepresentation;
        }

        return "c";
    }

    public boolean isHelixType() {
        return name().contains("HELIX");
    }

    public boolean isStrandType() {
        return equals(BRIDGE) || equals(EXTENDED);
    }

    public boolean isCoilType() {
        return !isHelixType() && !isStrandType();
    }
}
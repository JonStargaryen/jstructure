package de.bioforscher.jstructure.feature.sse;

/**
 * The enum of DSSP secondary structure elements.
 * Created by S on 04.10.2016.
 */
public enum SecondaryStructureElement {
    COIL("c"),
    BEND("S"),
    TURN("T"),
    PIHELIX("I"),
    THREE10HELIX("G"),
    BRIDGE("B"),
    EXTENDED("E"),
    ALPHA_HELIX("H");

    private String oneLetterRepresentation;

    SecondaryStructureElement(String oneLetterRepresentation) {
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

        return " ";
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
package de.bioforscher.jstructure.feature.sse;

/**
 * Created by S on 04.10.2016.
 */
public enum DSSPSecondaryStructureElement {
    COIL,
    BEND,
    TURN,
    PIHELIX,
    THREE10HELIX,
    BRIDGE,
    EXTENDED,
    ALPHA_HELIX;

    public boolean isHelixType() {
        return this.name().contains("HELIX");
    }

    public boolean isStrandType() {
        return this.equals(BRIDGE) || this.equals(EXTENDED);
    }
}
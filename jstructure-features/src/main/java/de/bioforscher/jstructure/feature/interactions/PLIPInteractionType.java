package de.bioforscher.jstructure.feature.interactions;

import java.util.Optional;
import java.util.stream.Stream;

/**
 * Enumerates all possible PLIP interaction types.
 * Created by bittrich on 2/15/17.
 */
public enum PLIPInteractionType {
    HALOGEN_BOND("halogen_bond", HalogenBond.class),
    HYDROGEN_BOND("hydrogen_bond", HydrogenBond.class),
    HYDROPHOBIC_INTERACTION("hydrophobic_interaction", HydrophobicInteraction.class),
    METAL_COMPLEX("metal_complex", MetalComplex.class),
    PI_CATION_INTERACTION("pi_cation_interaction", PiCationInteraction.class),
    PI_STACKING("pi_stack", PiStacking.class),
    SALT_BRIDGE("salt_bridge", SaltBridge.class),
    WATER_BRIDGE("water_bridge", WaterBridge.class);

    private String interactionTag;
    private Class<? extends PLIPInteraction> describingClass;

    PLIPInteractionType(String interactionTag, Class<? extends PLIPInteraction> describingClass) {
        this.interactionTag = interactionTag;
        this.describingClass = describingClass;
    }

    public static Optional<PLIPInteractionType> valueOf(Class<? extends PLIPInteraction> describingClass) {
        return Stream.of(values())
                .filter(interaction -> interaction.getDescribingClass().equals(describingClass))
                .findFirst();
    }

    public String getInteractionTag() {
        return interactionTag;
    }

    public Class<? extends PLIPInteraction> getDescribingClass() {
        return describingClass;
    }
}

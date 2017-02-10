package de.bioforscher.jstructure.parser.plip;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumerates all types of PLIP interactions.
 * Created by bittrich on 2/9/17.
 */
public enum PLIPInteractionType {
    HYDROGEN_BOND("HBond", "hydrogen_bond", "donoridx", "acceptoridx"),
    WATER_BRIDGE("WaterBr", "water_bridge", "donor_idx", "acceptor_idx"),
    METAL_COMPLEX("Metal", "metal_complex", "metal_idx", "target_idx"),
    SALT_BRIDGE("SaltBr", "salt_bridge", "lig_idx_list"),
    HYDROPHOBIC("Phobic", "hydrophobic_interaction", "ligcarbonidx", "protcarbonidx"),
    PI_STACKING("PiStack", "pi_stack", "lig_idx_list"),
    PI_CATION("PiCat", "pi_cation_interaction", "lig_idx_list"),
    HALOGEN_BOND("HaloBond", "halogen_bond", "don_idx", "acc_idx");

    /**
     * A more readable label.
     */
    private String shortName;
    /**
     * Which tag in the PLIP-xml contains this information?
     */
    private String interactionTag;
    /**
     * What to parse?
     */
    private List<String> atomTags;

    PLIPInteractionType(String shortName, String interactionTag, String... atomTags) {
        this.shortName = shortName;
        this.interactionTag = interactionTag;
        this.atomTags = Stream.of(atomTags).collect(Collectors.toList());
    }

    String getInteractionTag() {
        return interactionTag;
    }

    List<String> getAtomTags() {
        return atomTags;
    }

    String getShortName() {
        return shortName;
    }
}

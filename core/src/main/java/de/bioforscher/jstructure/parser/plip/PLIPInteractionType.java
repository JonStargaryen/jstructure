package de.bioforscher.jstructure.parser.plip;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumerates all types of PLIP interactions.
 * Created by bittrich on 2/9/17.
 */
public enum PLIPInteractionType {
    HYDROGEN_BOND("hydrogen_bond", "donoridx", "acceptoridx"),
    WATER_BRIDGE("water_bridge", "donor_idx", "acceptor_idx"),
    METAL_COMPLEX("metal_complex", "metal_idx", "target_idx"),
    SALT_BRIDGE("salt_bridge", "lig_idx_list"),
    HYDROPHOBIC("hydrophobic_interaction", "ligcarbonidx", "protcarbonidx"),
    PI_STACKING("pi_stack", "lig_idx_list"),
    PI_CATION("pi_cation_interaction", "lig_idx_list"),
    HALOGEN_BOND("halogen_bond", "don_idx", "acc_idx")
    ;

    private String interactionTag;
    private List<String> atomTags;

    PLIPInteractionType(String interactionTag, String... atomTags) {
        this.interactionTag = interactionTag;
        this.atomTags = Stream.of(atomTags).collect(Collectors.toList());
    }

    String getInteractionTag() {
        return interactionTag;
    }

    List<String> getAtomTags() {
        return atomTags;
    }
}

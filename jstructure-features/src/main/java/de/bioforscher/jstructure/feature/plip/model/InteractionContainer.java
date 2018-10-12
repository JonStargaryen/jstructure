package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class InteractionContainer {
    private final List<Interaction> interactions;
    private final List<Interaction> polymerInteractions;
    private final List<Interaction> ligandInteractions;

    private final List<HalogenBond> halogenBonds;
    private final List<HydrogenBond> hydrogenBonds;
    private final List<HydrophobicInteraction> hydrophobicInteractions;
    private final List<MetalComplex> metalComplexes;
    private final List<PiCationInteraction> piCationInteractions;
    private final List<PiStackingInteraction> piStackingInteractions;
    private final List<SaltBridge> saltBridges;
    private final List<WaterBridge> waterBridges;

    public InteractionContainer(List<Interaction> interactions) {
        this.interactions = interactions;
        this.polymerInteractions = extractInteractions(interactions, Interaction::isPolymerInteraction);
        this.ligandInteractions = extractInteractions(interactions, Interaction::isLigandInteraction);

        this.halogenBonds = extractAndCastInteractions(interactions, HalogenBond.class);
        this.hydrogenBonds = extractAndCastInteractions(interactions, HydrogenBond.class);
        this.hydrophobicInteractions = extractAndCastInteractions(interactions, HydrophobicInteraction.class);
        this.metalComplexes = extractAndCastInteractions(interactions, MetalComplex.class);
        this.piCationInteractions = extractAndCastInteractions(interactions, PiCationInteraction.class);
        this.piStackingInteractions = extractAndCastInteractions(interactions, PiStackingInteraction.class);
        this.saltBridges = extractAndCastInteractions(interactions, SaltBridge.class);
        this.waterBridges = extractAndCastInteractions(interactions, WaterBridge.class);
    }

    public List<Interaction> getInteractions() {
        return interactions;
    }

    public InteractionContainer getSubsetOfInteractions(Group group) {
        return new InteractionContainer(interactions.stream()
                .filter(interaction -> interaction.containsGroup(group))
                .collect(Collectors.toList()));
    }

    public List<Interaction> getPolymerInteractions() {
        return polymerInteractions;
    }

    public List<Interaction> getLigandInteractions() {
        return ligandInteractions;
    }

    public List<HalogenBond> getHalogenBonds() {
        return halogenBonds;
    }

    public List<HydrogenBond> getHydrogenBonds() {
        return hydrogenBonds;
    }

    public List<HydrophobicInteraction> getHydrophobicInteractions() {
        return hydrophobicInteractions;
    }

    public List<MetalComplex> getMetalComplexes() {
        return metalComplexes;
    }

    public List<PiCationInteraction> getPiCationInteractions() {
        return piCationInteractions;
    }

    public List<PiStackingInteraction> getPiStackingInteractions() {
        return piStackingInteractions;
    }

    public List<SaltBridge> getSaltBridges() {
        return saltBridges;
    }

    public List<WaterBridge> getWaterBridges() {
        return waterBridges;
    }

    private List<Interaction> extractInteractions(List<Interaction> interactions,
                                                  Predicate<Interaction> predicate) {
        return interactions.stream()
                .filter(predicate)
                .collect(Collectors.toList());
    }

    private <I extends AbstractInteraction> List<I> extractAndCastInteractions(List<Interaction> interactions,
                                                                       Class<I> targetClass) {
        return interactions.stream()
                .filter(targetClass::isInstance)
                .map(targetClass::cast)
                .collect(Collectors.toList());
    }
}

package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.Featureable;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A container for all {@link PLIPInteraction} instances present for a
 * given {@link Featureable}.
 * Created by bittrich on 2/16/17.
 */
public class PLIPInteractionContainer extends FeatureContainerEntry {
    private final List<HalogenBond> halogenBonds;
    private final List<HydrogenBond> hydrogenBonds;
    private final List<HydrophobicInteraction> hydrophobicInteractions;
    private final List<MetalComplex> metalComplexes;
    private final List<PiCationInteraction> piCationInteractions;
    private final List<PiStacking> piStackings;
    private final List<SaltBridge> saltBridges;
    private final List<WaterBridge> waterBridges;

    PLIPInteractionContainer(FeatureProvider featureProvider, List<PLIPInteraction> plipInteractions) {
        super(featureProvider);
        halogenBonds = filterInteractions(plipInteractions, HalogenBond.class);
        hydrogenBonds = filterInteractions(plipInteractions, HydrogenBond.class);
        hydrophobicInteractions = filterInteractions(plipInteractions, HydrophobicInteraction.class);
        metalComplexes = filterInteractions(plipInteractions, MetalComplex.class);
        piCationInteractions = filterInteractions(plipInteractions, PiCationInteraction.class);
        piStackings = filterInteractions(plipInteractions, PiStacking.class);
        saltBridges = filterInteractions(plipInteractions, SaltBridge.class);
        waterBridges = filterInteractions(plipInteractions, WaterBridge.class);
    }

    private <I> List<I> filterInteractions(List<PLIPInteraction> globalPlipInteractions, Class<I> content) {
        return globalPlipInteractions.stream()
                .filter(content::isInstance)
                .map(content::cast)
                .collect(Collectors.toList());
    }

    public List<PLIPInteraction> getInteractions() {
        return Stream.of(halogenBonds, hydrogenBonds, hydrophobicInteractions, metalComplexes, piCationInteractions,
                piStackings, saltBridges, waterBridges)
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
    }

    public List<PLIPInteraction> getBackboneInteractions() {
        return getInteractions().stream()
                .filter(PLIPInteraction::isBackboneInteraction)
                .collect(Collectors.toList());
    }

    public List<PLIPInteraction> getSideChainInteractions() {
        return getInteractions().stream()
                .filter(PLIPInteraction::isSideChainInteraction)
                .collect(Collectors.toList());
    }

    public List<PLIPInteraction> getMixedInteractions() {
        return getInteractions().stream()
                .filter(PLIPInteraction::isMixedInteraction)
                .collect(Collectors.toList());
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

    public List<PiStacking> getPiStackings() {
        return piStackings;
    }

    public List<SaltBridge> getSaltBridges() {
        return saltBridges;
    }

    public List<WaterBridge> getWaterBridges() {
        return waterBridges;
    }

    public PLIPInteractionContainer getInteractionsFor(Group group) {
        List<PLIPInteraction> filteredInteractions = getInteractions().stream()
                .filter(plipInteraction -> plipInteraction.getPartner1().equals(group) || plipInteraction.getPartner2().equals(group))
                .collect(Collectors.toList());
        return new PLIPInteractionContainer(getFeatureProvider(), filteredInteractions);
    }

    public boolean areInContact(Group group1, Group group2) {
        return getInteractions().stream()
                .anyMatch(plipInteraction -> (plipInteraction.getPartner1().equals(group1) && plipInteraction.getPartner2().equals(group2)) ||
                        (plipInteraction.getPartner2().equals(group1) && plipInteraction.getPartner1().equals(group2)));
    }
}

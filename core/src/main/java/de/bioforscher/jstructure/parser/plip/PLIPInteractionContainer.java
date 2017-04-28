package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.parser.plip.interaction.*;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A container for all {@link de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction} instances present for a
 * given {@link de.bioforscher.jstructure.model.feature.FeatureContainer}.
 * Created by bittrich on 2/16/17.
 */
public class PLIPInteractionContainer {
    private final List<HalogenBond> halogenBonds;
    private final List<HydrogenBond> hydrogenBonds;
    private final List<HydrophobicInteraction> hydrophobicInteractions;
    private final List<MetalComplex> metalComplexes;
    private final List<PiCationInteraction> piCationInteractions;
    private final List<PiStacking> piStackings;
    private final List<SaltBridge> saltBridges;
    private final List<WaterBridge> waterBridges;

    public PLIPInteractionContainer(List<PLIPInteraction> plipInteractions) {
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

    PLIPInteractionContainer(List<HalogenBond> halogenBonds, List<HydrogenBond> hydrogenBonds, List<HydrophobicInteraction> hydrophobicInteractions, List<MetalComplex> metalComplexes, List<PiCationInteraction> piCationInteractions, List<PiStacking> piStackings, List<SaltBridge> saltBridges, List<WaterBridge> waterBridges) {
        this.halogenBonds = halogenBonds;
        this.hydrogenBonds = hydrogenBonds;
        this.hydrophobicInteractions = hydrophobicInteractions;
        this.metalComplexes = metalComplexes;
        this.piCationInteractions = piCationInteractions;
        this.piStackings = piStackings;
        this.saltBridges = saltBridges;
        this.waterBridges = waterBridges;
    }

    public List<PLIPInteraction> getInteractions() {
        return Stream.of(halogenBonds, hydrogenBonds, hydrophobicInteractions, metalComplexes, piCationInteractions,
                piStackings, saltBridges, waterBridges)
                .flatMap(Collection::stream)
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
}

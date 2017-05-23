package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.Featureable;

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
    private final List<MetalComplex> metalComplexes;
    private final List<PiCationInteraction> piCationInteractions;
    private final List<PiStacking> piStackings;
    private final List<SaltBridge> saltBridges;
    private final List<WaterBridge> waterBridges;

    PLIPInteractionContainer(AbstractFeatureProvider featureProvider, List<PLIPInteraction> plipInteractions) {
        super(featureProvider);
        halogenBonds = filterInteractions(plipInteractions, HalogenBond.class);
        hydrogenBonds = filterInteractions(plipInteractions, HydrogenBond.class);
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
        return Stream.of(halogenBonds, hydrogenBonds, metalComplexes, piCationInteractions,
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

    //TODO checks for given residues and convenience functions to retrieve information
}

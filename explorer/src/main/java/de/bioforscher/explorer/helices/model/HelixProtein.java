package de.bioforscher.explorer.helices.model;

import de.bioforscher.explorer.membrane.model.ExplorerInteraction;
import de.bioforscher.explorer.membrane.model.ExplorerLigand;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 3/7/17.
 */
public class HelixProtein {
    private String name, title;
    private List<HelixChain> chains;
    /* interactions */
    private List<ExplorerInteraction> halogenBonds;
    private List<ExplorerInteraction> hydrogenBonds;
    private List<ExplorerInteraction> metalComplexes;
    private List<ExplorerInteraction> piCationInteractions;
    private List<ExplorerInteraction> piStackings;
    private List<ExplorerInteraction> saltBridges;
    private List<ExplorerInteraction> waterBridges;
    /* ligands */
    private List<ExplorerLigand> ligands;

    public HelixProtein() {

    }

    private static HelixProtein createHelixProteinFromScaffold(Protein protein) {
        HelixProtein explorerProtein = new HelixProtein();

        explorerProtein.name = protein.getName();
        explorerProtein.title = protein.getTitle();
        explorerProtein.chains = protein.chains()
                .map(HelixChain::new)
                .collect(Collectors.toList());

        explorerProtein.ligands = Selection.on(protein)
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());

        return explorerProtein;
    }

    public static HelixProtein createForHelixExplorer(String pdbId) {
        Protein protein = ProteinParser.parseProteinById(pdbId);
        AbstractFeatureProvider plipAnnotator = FeatureProviderRegistry.resolveAnnotator(PLIPAnnotator.PLIP_INTERACTIONS);
        plipAnnotator.process(protein);

        HelixProtein helixProtein = createHelixProteinFromScaffold(protein);
        PLIPInteractionContainer interactions = protein.getFeature(PLIPInteractionContainer.class, PLIPAnnotator.PLIP_INTERACTIONS);
        helixProtein.halogenBonds = interactions.getHalogenBonds().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        helixProtein.hydrogenBonds = interactions.getHydrogenBonds().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        helixProtein.metalComplexes = interactions.getMetalComplexes().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        helixProtein.piCationInteractions = interactions.getPiCationInteractions().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        helixProtein.piStackings = interactions.getPiStackings().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        helixProtein.saltBridges = interactions.getSaltBridges().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        helixProtein.waterBridges = interactions.getWaterBridges().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());

        return helixProtein;
    }

    public String getName() {
        return name;
    }

    public String getTitle() {
        return title;
    }

    public List<HelixChain> getChains() {
        return chains;
    }

    public List<ExplorerInteraction> getHalogenBonds() {
        return halogenBonds;
    }

    public List<ExplorerInteraction> getHydrogenBonds() {
        return hydrogenBonds;
    }

    public List<ExplorerInteraction> getMetalComplexes() {
        return metalComplexes;
    }

    public List<ExplorerInteraction> getPiCationInteractions() {
        return piCationInteractions;
    }

    public List<ExplorerInteraction> getPiStackings() {
        return piStackings;
    }

    public List<ExplorerInteraction> getSaltBridges() {
        return saltBridges;
    }

    public List<ExplorerInteraction> getWaterBridges() {
        return waterBridges;
    }

    public List<ExplorerLigand> getLigands() {
        return ligands;
    }
}

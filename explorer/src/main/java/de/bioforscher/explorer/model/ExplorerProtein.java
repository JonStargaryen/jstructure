package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderHelix;
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The reduced representation of {@link Protein} objects
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerProtein {
    private String name, uniprot, title;
    private List<ExplorerChain> chains;
    /* sequence motifs */
    private List<ExplorerMotif> sequenceMotifs;
    /* topology information */
    private List<double[]> membrane;
    /* interactions */
    private List<ExplorerInteraction> halogenBonds;
    private List<ExplorerInteraction> hydrogenBonds;
    private List<ExplorerInteraction> hydrophobicInteractions;
    private List<ExplorerInteraction> metalComplexes;
    private List<ExplorerInteraction> piCationInteractions;
    private List<ExplorerInteraction> piStackings;
    private List<ExplorerInteraction> saltBridges;
    private List<ExplorerInteraction> waterBridges;
    /* ligands */
    private List<ExplorerLigand> ligands;
    /* helices & kinks */
    private List<ExplorerHelix> helices;
    private List<ExplorerKink> kinks;
    /* UniProt annotations are chain-specific and stored there */

    public ExplorerProtein() {

    }

    public ExplorerProtein(Protein protein) {
        //TODO standardize to helix explorer
        this.name = protein.getName();
        //TODO handle how uniprot ids are linked to chains
        this.uniprot = protein.chains()
                .map(chain -> chain.getFeature(String.class, UniProtAnnotator.UNIPROT_ID))
                .distinct()
                .collect(Collectors.joining(", "));
        this.title = protein.getTitle();
        this.chains = protein.chains()
                .map(ExplorerChain::new)
                .collect(Collectors.toList());

        this.sequenceMotifs = protein.getFeatureAsList(SequenceMotif.class, SequenceMotifAnnotator.SEQUENCE_MOTIF).stream()
                .map(ExplorerMotif::new)
                .collect(Collectors.toList());

        this.membrane = protein.getFeature(Membrane.class, ANVIL.MEMBRANE).getMembraneAtoms();

        PLIPInteractionContainer interactions = protein.getFeature(PLIPInteractionContainer.class, PLIPAnnotator.PLIP_INTERACTIONS);
        this.halogenBonds = interactions.getHalogenBonds().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        this.hydrogenBonds = interactions.getHydrogenBonds().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        this.hydrophobicInteractions = interactions.getHydrophobicInteractions().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        this.metalComplexes = interactions.getMetalComplexes().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        this.piCationInteractions = interactions.getPiCationInteractions().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                            .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        this.piStackings = interactions.getPiStackings().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        this.saltBridges = interactions.getSaltBridges().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        this.waterBridges = interactions.getWaterBridges().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());

        this.ligands = Selection.on(protein)
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());

        this.helices = protein.chains()
                .flatMap(ExplorerHelix::extract)
                .collect(Collectors.toList());

        this.kinks = protein.getFeatureAsList(KinkFinderHelix.class, KinkFinderParser.KINK_FINDER_ANNOTATION).stream()
                .map(ExplorerKink::new)
                .collect(Collectors.toList());
    }

    private static ExplorerProtein createExplorerProteinFromScaffold(Protein protein) {
        ExplorerProtein explorerProtein = new ExplorerProtein();

        explorerProtein.name = protein.getName();
        explorerProtein.title = protein.getTitle();
        explorerProtein.chains = protein.chains()
                .map(ExplorerChain::new)
                .collect(Collectors.toList());

        explorerProtein.ligands = Selection.on(protein)
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());

        return explorerProtein;
    }

    public static ExplorerProtein createForHelixExplorer(String pdbId) {
        Protein protein = ProteinParser.parseProteinById(pdbId);
        AbstractFeatureProvider plipAnnotator = FeatureProviderRegistry.resolveAnnotator(PLIPAnnotator.PLIP_INTERACTIONS);
        plipAnnotator.process(protein);

        ExplorerProtein explorerProtein = createExplorerProteinFromScaffold(protein);
        PLIPInteractionContainer interactions = protein.getFeature(PLIPInteractionContainer.class, PLIPAnnotator.PLIP_INTERACTIONS);
        explorerProtein.halogenBonds = interactions.getHalogenBonds().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        explorerProtein.hydrogenBonds = interactions.getHydrogenBonds().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        explorerProtein.metalComplexes = interactions.getMetalComplexes().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());
        explorerProtein.piCationInteractions = interactions.getPiCationInteractions().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        explorerProtein.piStackings = interactions.getPiStackings().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        explorerProtein.saltBridges = interactions.getSaltBridges().stream()
                .filter(interaction -> !interaction.getAtoms1().isEmpty())
                .flatMap(interaction -> Combinatorics.cartesianProductOf(interaction.getAtoms1(), interaction.getAtoms2())
                        .map(ExplorerInteraction::new))
                .collect(Collectors.toList());
        explorerProtein.waterBridges = interactions.getWaterBridges().stream()
                .map(ExplorerInteraction::new)
                .collect(Collectors.toList());

        return explorerProtein;
    }

    public String getName() {
        return name;
    }

    public String getUniprot() {
        return uniprot;
    }

    public String getTitle() {
        return title;
    }

    public List<ExplorerChain> getChains() {
        return chains;
    }

    public List<ExplorerMotif> getSequenceMotifs() {
        return sequenceMotifs;
    }

    public List<double[]> getMembrane() {
        return membrane;
    }

    public List<ExplorerInteraction> getHalogenBonds() {
        return halogenBonds;
    }

    public List<ExplorerInteraction> getHydrogenBonds() {
        return hydrogenBonds;
    }

    public List<ExplorerInteraction> getHydrophobicInteractions() {
        return hydrophobicInteractions;
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

    public List<ExplorerHelix> getHelices() {
        return helices;
    }

    public List<ExplorerKink> getKinks() {
        return kinks;
    }
}

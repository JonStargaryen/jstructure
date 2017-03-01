package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderHelix;
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The reduced representation of {@link Protein} objects
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerProtein {
    private String name, title;
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

    public ExplorerProtein() {

    }

    public ExplorerProtein(Protein protein) {
        this.name = protein.getName();
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

    public String getName() {
        return name;
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

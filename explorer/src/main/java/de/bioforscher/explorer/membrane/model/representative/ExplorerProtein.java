package de.bioforscher.explorer.membrane.model.representative;

import de.bioforscher.explorer.membrane.model.homologous.HomologousProteinContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderHelix;
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderParser;
import de.bioforscher.jstructure.parser.opm.OPMDatabaseQuery;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;

import java.util.ArrayList;
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
    /* homologous protein ids for this structure */
    private HomologousProteinContainer homologous;
    private boolean noPlipData;
    /* UniProt annotations are chain-specific and stored there */

    public ExplorerProtein() {

    }

    public ExplorerProtein(Protein protein) {
        this.name = protein.getName();
        this.title = protein.getTitle();
        this.chains = protein.chains()
                .map(ExplorerChain::new)
                .collect(Collectors.toList());

        try {
            this.sequenceMotifs = protein.getFeatureAsList(SequenceMotif.class, SequenceMotifAnnotator.SEQUENCE_MOTIF).stream()
                    .map(ExplorerMotif::new)
                    .collect(Collectors.toList());
        } catch (NullPointerException e) {
            this.sequenceMotifs = new ArrayList<>();
        }

        try {
            this.membrane = protein.getFeature(Membrane.class, ANVIL.MEMBRANE).getMembraneAtoms();
        } catch (NullPointerException e) {
            this.membrane = new ArrayList<>();
        }

        try {
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
        } catch (NullPointerException e) {
            this.halogenBonds = new ArrayList<>();
            this.hydrogenBonds = new ArrayList<>();
            this.hydrophobicInteractions = new ArrayList<>();
            this.metalComplexes = new ArrayList<>();
            this.piCationInteractions = new ArrayList<>();
            this.piStackings = new ArrayList<>();
            this.saltBridges = new ArrayList<>();
            this.waterBridges = new ArrayList<>();
            this.noPlipData = true;
        }

        this.ligands = Selection.on(protein)
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());

        try {
            this.helices = protein.chains()
                    .flatMap(ExplorerHelix::extract)
                    .collect(Collectors.toList());
        } catch (NullPointerException e) {
            this.helices = new ArrayList<>();
        }

        try {
            this.kinks = protein.getFeatureAsList(KinkFinderHelix.class, KinkFinderParser.KINK_FINDER_ANNOTATION).stream()
                    .map(ExplorerKink::new)
                    .collect(Collectors.toList());
        } catch (NullPointerException e) {
            this.kinks = new ArrayList<>();
        }

        this.homologous = new HomologousProteinContainer(protein.getFeatureAsList(String.class, OPMDatabaseQuery.HOMOLOGOUS_PROTEINS));
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

    public HomologousProteinContainer getHomologous() {
        return homologous;
    }

    public boolean isNoPlipData() {
        return noPlipData;
    }
}

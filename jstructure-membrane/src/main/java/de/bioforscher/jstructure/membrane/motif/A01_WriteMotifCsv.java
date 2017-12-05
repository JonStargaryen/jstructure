package de.bioforscher.jstructure.membrane.motif;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.stream.Collectors;

/**
 * Evaluate all sequence motifs in the data set.
 */
public class A01_WriteMotifCsv {
    private static final Path DIRECTORY = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;
    private static final SequenceMotifAnnotator SEQUENCE_MOTIF_ANNOTATOR = new SequenceMotifAnnotator();

    public static void main(String[] args) {
        // general stuff
        String header = "pdbId;chainId;startResId;motif;sequence;topology;" +
                // DynaMine stuff
                "rigidity;" +
                // evolutionary information
                "evol;" +
                // backbone interactions
                "bb;" +
                // non-local interactions - at least 6 residues apart
                "nl_halogen;nl_hydrogen;nl_hydrophobic;nl_metal;nl_picat;nl_pistack;nl_salt;nl_water;nl_total;" +
                // helix-helix interactions - different helix, >15 residues
                "hh_halogen;hh_hydrogen;hh_hydrophobic;hh_metal;hh_picat;hh_pistack;hh_salt;hh_water;hh_total" + System.lineSeparator();

        String output = MembraneConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A01_WriteMotifCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(DIRECTORY.resolve("results").resolve("motifs.csv"),
                header + output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            Optional<MembraneConstants.WrappedChain> wrappedChainOptional = MembraneConstants.handleLine(line, DIRECTORY);

            if(!wrappedChainOptional.isPresent()) {
                return Optional.empty();
            }

            MembraneConstants.WrappedChain wrappedChain = wrappedChainOptional.get();
            Chain chain = wrappedChain.getChain();
            String pdbId = wrappedChain.getPdbId();
            String chainId = wrappedChain.getChainId();
            SEQUENCE_MOTIF_ANNOTATOR.process(chain.getParentStructure());
            boolean evolScoresSane = wrappedChain.isEvolScoresSane();

            String output = chain.getFeature(SequenceMotifContainer.class)
                    .getSequenceMotifs()
                    .stream()
                    .map(sequenceMotif -> handleSequenceMotif(pdbId,
                            chainId,
                            sequenceMotif,
                            evolScoresSane))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator()));

            return Optional.of(output);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleSequenceMotif(String pdbId,
                                                    String chainId,
                                                    SequenceMotif sequenceMotif,
                                                    boolean evolScoresSane) {
        try {
            String topology = determineTopology(sequenceMotif);
            boolean tm = topology.equals("I");

            List<AminoAcid> aminoAcids = sequenceMotif.getAminoAcids();

            OptionalDouble backboneRigidity = aminoAcids.stream()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeatureOptional(BackboneRigidity.class))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .mapToDouble(BackboneRigidity::getBackboneRigidity)
                    .average();
            OptionalDouble evolScore = aminoAcids.stream()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeatureOptional(EvolutionaryInformation.class))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .mapToDouble(EvolutionaryInformation::getInformation)
                    .average();

            PLIPInteractionContainer interactions = new PLIPInteractionContainer(null,
                    aminoAcids.stream()
                            .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                            .map(PLIPInteractionContainer::getInteractions)
                            .flatMap(Collection::stream)
                            .collect(Collectors.toList()));

            // evaluate interaction types
            PLIPInteractionContainer nonLocalInteractions = new PLIPInteractionContainer(null,
                    interactions.getInteractions()
                            .stream()
                            // interactions have to be non-local
                            .filter(interaction -> Math.abs(interaction.getPartner1().getResidueIdentifier().getResidueNumber() - interaction.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                            .collect(Collectors.toList()));
            // extract all interaction with helices - either other TM-helices or other helices, >15 residues
            PLIPInteractionContainer helixHelixInteractions = new PLIPInteractionContainer(null,
                    nonLocalInteractions.getInteractions()
                            .stream()
                            .filter(interaction -> tm ? MembraneConstants.isTmHelixInteraction(interaction) : MembraneConstants.isNtmHelixInteraction(interaction))
                            .collect(Collectors.toList()));

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    sequenceMotif.getAminoAcids().get(0).getResidueIdentifier().getResidueNumber() + ";" +
                    sequenceMotif.getMotifDefinition().name() + ";" +
                    aminoAcids.stream().map(AminoAcid::getOneLetterCode).collect(Collectors.joining()) + ";" +
                    topology + ";" +

                    (backboneRigidity.isPresent() ? StandardFormat.format(backboneRigidity.getAsDouble()) : "NA") + ";" +

                    (evolScoresSane && evolScore.isPresent() ? StandardFormat.format(evolScore.getAsDouble()) : "NA") + ";" +

                    interactions.getBackboneInteractions().size() + ";" +

                    nonLocalInteractions.getHalogenBonds().size() + ";" +
                    nonLocalInteractions.getHydrogenBonds().size() + ";" +
                    nonLocalInteractions.getHydrophobicInteractions().size() + ";" +
                    nonLocalInteractions.getMetalComplexes().size() + ";" +
                    nonLocalInteractions.getPiCationInteractions().size() + ";" +
                    nonLocalInteractions.getPiStackings().size() + ";" +
                    nonLocalInteractions.getSaltBridges().size() + ";" +
                    nonLocalInteractions.getWaterBridges().size() + ";" +
                    nonLocalInteractions.getInteractions().size() + ";" +

                    helixHelixInteractions.getHalogenBonds().size() + ";" +
                    helixHelixInteractions.getHydrogenBonds().size() + ";" +
                    helixHelixInteractions.getHydrophobicInteractions().size() + ";" +
                    helixHelixInteractions.getMetalComplexes().size() + ";" +
                    helixHelixInteractions.getPiCationInteractions().size() + ";" +
                    helixHelixInteractions.getPiStackings().size() + ";" +
                    helixHelixInteractions.getSaltBridges().size() + ";" +
                    helixHelixInteractions.getWaterBridges().size() + ";" +
                    helixHelixInteractions.getInteractions().size());
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static String determineTopology(SequenceMotif sequenceMotif) {
        if(sequenceMotif.getAminoAcids()
            .stream()
            .allMatch(aminoAcid -> aminoAcid.getFeature(Topology.class).isTransmembrane())) {
            return "I";
        }
        if(sequenceMotif.getAminoAcids()
                .stream()
                .noneMatch(aminoAcid -> aminoAcid.getFeature(Topology.class).isTransmembrane())) {
            return "o";
        }
        return "t";
    }
}

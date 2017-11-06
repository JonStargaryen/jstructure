package de.bioforscher.jstructure.membrane.aminoacid;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.Kink;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Evaluates all amino acids in the alpha-dataset. Indirectly covers interactions.
 */
public class A01_WriteAminoAcidContactCsv {
    private static final Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;

    public static void main(String[] args) {
                        // general stuff
        String header = "pdbId;chainId;resId;aa;sequence;sse3;sse9;sseSize;topology;" +
                // DynaMine stuff
                "rigidity;" +
                // KinkFinder stuff
                "kink;significantKink;kinkAngle;" +
                // evolutionary information
                "evol;" +
                // backbone interactions
                "bb;" +
                // non-local interactions - at least 6 residues apart
                "nl_halogen;nl_hydrogen;nl_hydrophobic;nl_metal;nl_picat;nl_pistack;nl_salt;nl_water;nl_total;" +
                // helix-helix interactions - different helix, >15 residues
                "hh_halogen;hh_hydrogen;hh_hydrophobic;hh_metal;hh_picat;hh_pistack;hh_salt;hh_water;hh_total" + System.lineSeparator();

        String output = MembraneConstants.lines(directory.resolve("ids.list"))
                .map(A01_WriteAminoAcidContactCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(directory.resolve("results").resolve("aminoacids.csv"),
                header + output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            Optional<MembraneConstants.WrappedChain> wrappedChainOptional = MembraneConstants.handleLine(line, directory);

            if(!wrappedChainOptional.isPresent()) {
                return Optional.empty();
            }

            MembraneConstants.WrappedChain wrappedChain = wrappedChainOptional.get();
            Chain chain = wrappedChain.getChain();
            String pdbId = wrappedChain.getPdbId();
            String chainId = wrappedChain.getChainId();
            boolean evolScoresSane = wrappedChain.isEvolScoresSane();
            boolean kinkFinderDataSane = wrappedChain.isKinkFinderDataSane();

            String output = chain.aminoAcids()
                    .map(aminoAcid -> handleAminoAcid(pdbId,
                            chainId,
                            aminoAcid,
                            evolScoresSane,
                            kinkFinderDataSane))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator()));

            return Optional.of(output);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleAminoAcid(String pdbId,
                                                    String chainId,
                                                    AminoAcid aminoAcid,
                                                    boolean evolScoresSane,
                                                    boolean kinkFinderDataSane) {
        try {
            boolean tm = aminoAcid.getFeature(Topology.class).isTransmembrane();
            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement =
                    aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid);

            PLIPInteractionContainer interactions = aminoAcid.getFeature(PLIPInteractionContainer.class);

            // evaluate interaction types
            PLIPInteractionContainer nonLocalInteractions = new PLIPInteractionContainer(null,
                    interactions
                            .getInteractions()
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

            GenericSecondaryStructure secondaryStructure = aminoAcid.getFeature(GenericSecondaryStructure.class);
            Optional<Kink> kink = aminoAcid.getFeatureContainer().getFeatureOptional(Kink.class);

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    aminoAcid.getResidueIdentifier() + ";" +
                    aminoAcid.getOneLetterCode() + ";" +
                    MembraneConstants.surroundingSequence(aminoAcid) + ";" +
                    secondaryStructure.getSecondaryStructure().getReducedRepresentation() + ";" +
                    secondaryStructure.getSecondaryStructure().getOneLetterRepresentation() + ";" +
                    surroundingSecondaryStructureElement.getSize() + ";" +
                    (tm ? "I" : "o") + ";" +

                    aminoAcid.getFeatureContainer()
                            .getFeatureOptional(BackboneRigidity.class)
                            .map(BackboneRigidity::getBackboneRigidity)
                            .map(StandardFormat::format)
                            .orElse("NA") + ";" +

                    (kinkFinderDataSane ? (kink.isPresent() ? "kinked" : "normal") : "NA") + ";" +
                    (kinkFinderDataSane ? (kink.map(kink1 -> (kink1.isSignificantKink() ? "signficiant" : "insignificant")).orElse("NA")) : "NA") + ";" +
                    (kinkFinderDataSane ? (kink.isPresent() ? StandardFormat.format(kink.get().getAngle()) : "NA") : "NA") + ";" +

                    (evolScoresSane ? StandardFormat.format(aminoAcid.getFeature(EvolutionaryInformation.class).getInformation()) : "NA") + ";" +

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
}

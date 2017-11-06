package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.membrane.Kink;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Evaluates all interactions in the alpha-dataset. Indirectly covers amino acids.
 */
public class A01_WriteContactCsv {
    private static final Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;

    public static void main(String[] args) {
                        // general stuff
        String header = "pdbId;chainId;resId1;resId2;aa1;aa2;sequence1;sequence2;sse31;sse91;sse32;sse92;sseSize1;sseSize2;" +
                // interaction specific stuff
                "type;distance;" +
                // DynaMine stuff
                "rigidity1;rigidity2;" +
                // KinkFinder stuff
                "kink1;significantKink1;kinkAngle1;kink2;significantKink2;kinkAngle2;" +
                // evolutionary information
                "evol1;evol2;" +
                // backbone interactions
                "bb1;bb2;" +
                // non-local interactions - at least 6 residues apart
                "nl_halogen1;nl_hydrogen1;nl_hydrophobic1;nl_metal1;nl_picat1;nl_pistack1;nl_salt1;nl_water1;nl_total1;" +
                // helix-helix interactions - different helix, >15 residues
                "nl_halogen2;nl_hydrogen2;nl_hydrophobic2;nl_metal2;nl_picat2;nl_pistack2;nl_salt2;nl_water2;nl_total2;" + System.lineSeparator();

        String output = MembraneConstants.lines(directory.resolve("ids.list"))
                .map(A01_WriteContactCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(directory.resolve("results").resolve("contacts.csv"),
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

            String output = chain.getFeature(PLIPInteractionContainer.class)
                    .getInteractions()
                    .stream()
                    .filter(MembraneConstants::isTmHelixInteraction)
                    .map(interaction -> handleInteraction(pdbId,
                            chainId,
                            interaction,
                            evolScoresSane,
                            kinkFinderDataSane))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator()));

            if(output.isEmpty()) {
                return Optional.empty();
            }

            return Optional.of(output);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleInteraction(String pdbId,
                                                      String chainId,
                                                      PLIPInteraction interaction,
                                                      boolean evolScoresSane,
                                                      boolean kinkFinderDataSane) {
        try {
            AminoAcid aminoAcid1 = (AminoAcid) interaction.getPartner1();
            AminoAcid aminoAcid2 = (AminoAcid) interaction.getPartner2();

            PLIPInteractionContainer interactions1 = aminoAcid1.getFeature(PLIPInteractionContainer.class);
            PLIPInteractionContainer interactions2 = aminoAcid2.getFeature(PLIPInteractionContainer.class);

            // evaluate interaction types
            PLIPInteractionContainer nonLocalInteractions1 = new PLIPInteractionContainer(null,
                    interactions1
                            .getInteractions()
                            .stream()
                            // interactions have to be non-local
                            .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                            .collect(Collectors.toList()));
            PLIPInteractionContainer nonLocalInteractions2 = new PLIPInteractionContainer(null,
                    interactions2
                            .getInteractions()
                            .stream()
                            // interactions have to be non-local
                            .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                            .collect(Collectors.toList()));

            GenericSecondaryStructure secondaryStructure1 = aminoAcid1.getFeature(GenericSecondaryStructure.class);
            Optional<Kink> kink1 = aminoAcid1.getFeatureContainer().getFeatureOptional(Kink.class);
            GenericSecondaryStructure secondaryStructure2 = aminoAcid2.getFeature(GenericSecondaryStructure.class);
            Optional<Kink> kink2 = aminoAcid2.getFeatureContainer().getFeatureOptional(Kink.class);
            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement1 =
                    aminoAcid1.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid1);
            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement2 =
                    aminoAcid2.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid2);

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    aminoAcid1.getResidueIdentifier() + ";" +
                    aminoAcid2.getResidueIdentifier() + ";" +
                    aminoAcid1.getOneLetterCode() + ";" +
                    aminoAcid2.getOneLetterCode() + ";" +
                    MembraneConstants.surroundingSequence(aminoAcid1) + ";" +
                    MembraneConstants.surroundingSequence(aminoAcid2) + ";" +
                    secondaryStructure1.getSecondaryStructure().getReducedRepresentation() + ";" +
                    secondaryStructure1.getSecondaryStructure().getOneLetterRepresentation() + ";" +
                    secondaryStructure2.getSecondaryStructure().getReducedRepresentation() + ";" +
                    secondaryStructure2.getSecondaryStructure().getOneLetterRepresentation() + ";" +
                    surroundingSecondaryStructureElement1.getSize() + ";" +
                    surroundingSecondaryStructureElement2.getSize() + ";" +

                    interaction.getClass().getSimpleName() + ";" +
                    StandardFormat.format(aminoAcid1.getCa().calculate().distance(aminoAcid2.getCa())) + ";" +

                    aminoAcid1.getFeatureContainer()
                            .getFeatureOptional(BackboneRigidity.class)
                            .map(BackboneRigidity::getBackboneRigidity)
                            .map(StandardFormat::format)
                            .orElse("NA") + ";" +
                    aminoAcid2.getFeatureContainer()
                            .getFeatureOptional(BackboneRigidity.class)
                            .map(BackboneRigidity::getBackboneRigidity)
                            .map(StandardFormat::format)
                            .orElse("NA") + ";" +

                    (kinkFinderDataSane ? (kink1.isPresent() ? "kinked" : "normal") : "NA") + ";" +
                    (kinkFinderDataSane ? (kink1.map(k -> (k.isSignificantKink() ? "signficiant" : "insignificant")).orElse("NA")) : "NA") + ";" +
                    (kinkFinderDataSane ? (kink1.isPresent() ? StandardFormat.format(kink1.get().getAngle()) : "NA") : "NA") + ";" +
                    (kinkFinderDataSane ? (kink2.isPresent() ? "kinked" : "normal") : "NA") + ";" +
                    (kinkFinderDataSane ? (kink2.map(k -> (k.isSignificantKink() ? "signficiant" : "insignificant")).orElse("NA")) : "NA") + ";" +
                    (kinkFinderDataSane ? (kink2.isPresent() ? StandardFormat.format(kink2.get().getAngle()) : "NA") : "NA") + ";" +

                    (evolScoresSane ? StandardFormat.format(aminoAcid1.getFeature(EvolutionaryInformation.class).getInformation()) : "NA") + ";" +
                    (evolScoresSane ? StandardFormat.format(aminoAcid2.getFeature(EvolutionaryInformation.class).getInformation()) : "NA") + ";" +

                    interactions1.getBackboneInteractions().size() + ";" +
                    interactions2.getBackboneInteractions().size() + ";" +

                    nonLocalInteractions1.getHalogenBonds().size() + ";" +
                    nonLocalInteractions1.getHydrogenBonds().size() + ";" +
                    nonLocalInteractions1.getHydrophobicInteractions().size() + ";" +
                    nonLocalInteractions1.getMetalComplexes().size() + ";" +
                    nonLocalInteractions1.getPiCationInteractions().size() + ";" +
                    nonLocalInteractions1.getPiStackings().size() + ";" +
                    nonLocalInteractions1.getSaltBridges().size() + ";" +
                    nonLocalInteractions1.getWaterBridges().size() + ";" +
                    nonLocalInteractions1.getInteractions().size() + ";" +

                    nonLocalInteractions2.getHalogenBonds().size() + ";" +
                    nonLocalInteractions2.getHydrogenBonds().size() + ";" +
                    nonLocalInteractions2.getHydrophobicInteractions().size() + ";" +
                    nonLocalInteractions2.getMetalComplexes().size() + ";" +
                    nonLocalInteractions2.getPiCationInteractions().size() + ";" +
                    nonLocalInteractions2.getPiStackings().size() + ";" +
                    nonLocalInteractions2.getSaltBridges().size() + ";" +
                    nonLocalInteractions2.getWaterBridges().size() + ";" +
                    nonLocalInteractions2.getInteractions().size());
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}

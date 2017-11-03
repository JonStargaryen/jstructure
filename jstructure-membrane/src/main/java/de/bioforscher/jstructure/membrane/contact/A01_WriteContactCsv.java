package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.Kink;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Evaluates all interactions in the alpha-dataset. Indirectly covers amino acids.
 */
public class A01_WriteContactCsv {
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE =
            new DynaMineBridge();
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
            System.out.println("annotating " + line);
            String pdbId = line.split("_")[0];
            String chainId = line.split("_")[1];

            Structure structure = StructureParser.source(directory.resolve("pdb").resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();


            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                    Jsoup.parse(MembraneConstants.lines(directory.resolve("opm").resolve(pdbId + ".opm"))
                            .collect(Collectors.joining(System.lineSeparator()))));

            // skip chains with no or a single tm region
            int tmRegions = 0;
            boolean inTmRegion = false;
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            for (int i = 0; i < aminoAcids.size(); i++) {
                boolean tm = aminoAcids.get(i).getFeature(Topology.class).isTransmembrane();
                if(i == 0) {
                    inTmRegion = tm;
                    continue;
                }

                if(inTmRegion && !tm) {
                    tmRegions++;
                }

                inTmRegion = tm;
            }
            if(inTmRegion) {
                tmRegions++;
            }
            if(tmRegions < 2) {
                System.out.println("skipping structure with " + tmRegions + " TM-regions");
            }

            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                    Jsoup.parse(MembraneConstants.lines(directory.resolve("plip").resolve(line + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            DYNA_MINE_BRIDGE.process(chain,
                    MembraneConstants.lines(directory.resolve("dynamine").resolve(line + "_backbone.pred"))
                            .collect(Collectors.joining(System.lineSeparator())));
            List<Double> evolutionaryScores = MembraneConstants.lines(directory.resolve("evol").resolve(line + ".evol"))
                    .mapToDouble(Double::valueOf)
                    .boxed()
                    .collect(Collectors.toList());
            boolean evolScoresSane = aminoAcids.size() == evolutionaryScores.size();
            if(evolScoresSane) {
                for(int i = 0; i < aminoAcids.size(); i++) {
                    aminoAcids.get(i).getFeatureContainer().addFeature(new EvolutionaryInformation(null, null, evolutionaryScores.get(i)));
                }
            }

            // assign kink data
            boolean kinkException = false;
            try {
                MembraneConstants.lines(directory.resolve("kinks").resolve(pdbId + ".kinks"))
                        .filter(kinkLine -> !kinkLine.startsWith("pdb_code"))
                        .filter(kinkLine -> kinkLine.substring(4, 5).equals(chainId))
                        .forEach(kinkLine -> {
                            String[] split = kinkLine.split(",");
                            int kinkPosition = Integer.valueOf(split[3]);
                            double kinkAngle = Double.valueOf(split[6]);
                            boolean significantKink = kinkAngle > 15;
                            chain.select()
                                    .aminoAcids()
                                    .residueNumber(kinkPosition)
                                    .asAminoAcid()
                                    .getFeatureContainer()
                                    .addFeature(new Kink(kinkAngle, significantKink));
                        });
            } catch (Exception e) {
                e.printStackTrace();
                kinkException = true;
            }
            boolean kinkFinderDataSane = !kinkException;

            double min = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .min()
                    .orElse(0);
            double max = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .max()
                    .orElse(1);

            String output = chain.getFeature(PLIPInteractionContainer.class)
                    .getInteractions()
                    .stream()
                    .filter(MembraneConstants::isTmHelixInteraction)
                    .map(interaction -> handleInteraction(pdbId,
                            chainId,
                            interaction,
                            min,
                            max,
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

    private static Optional<String> handleInteraction(String pdbId,
                                                      String chainId,
                                                      PLIPInteraction interaction,
                                                      double min,
                                                      double max,
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
                            .map(value -> MembraneConstants.minMaxNormalize(value, min, max))
                            .map(StandardFormat::format)
                            .orElse("NA") + ";" +
                    aminoAcid2.getFeatureContainer()
                            .getFeatureOptional(BackboneRigidity.class)
                            .map(BackboneRigidity::getBackboneRigidity)
                            .map(value -> MembraneConstants.minMaxNormalize(value, min, max))
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

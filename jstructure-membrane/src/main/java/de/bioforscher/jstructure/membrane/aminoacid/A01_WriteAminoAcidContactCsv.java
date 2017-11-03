package de.bioforscher.jstructure.membrane.aminoacid;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
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
 * Evaluates all amino acids in the alpha-dataset. Indirectly covers interactions.
 */
public class A01_WriteAminoAcidContactCsv {
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE =
            new DynaMineBridge();
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

            String output = chain.aminoAcids()
                    .map(aminoAcid -> handleAminoAcid(pdbId,
                            chainId,
                            aminoAcid,
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

    private static Optional<String> handleAminoAcid(String pdbId,
                                                    String chainId,
                                                    AminoAcid aminoAcid,
                                                    double min,
                                                    double max,
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
                            .map(value -> MembraneConstants.minMaxNormalize(value, min, max))
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

package de.bioforscher.jstructure.efr.collect.obsolete;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interaction.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.graph.ResidueGraph;
import de.bioforscher.jstructure.graph.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.stream.Collectors;

public class S02_WriteMisfoldedNativeCsv {
    private static final Logger logger = LoggerFactory.getLogger(S02_WriteMisfoldedNativeCsv.class);
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.DATA_DIRECTORY.resolve("native").resolve("wozniak.list"))
                .map(S02_WriteMisfoldedNativeCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "exposed," +
                                "plip_hydrogen,plip_hydrophobic,plip_pistack,plip_cation,plip_salt,plip_metal,plip_halogen,plip_water,plip_bb,plip_total," +
                                "plip_l_hydrogen,plip_l_hydrophobic,plip_l_pistack,plip_l_cation,plip_l_salt,plip_l_metal,plip_halogen,plip_water,plip_l_bb,plip_l_total," +
                                "plip_nl_hydrogen,plip_nl_hydrophobic,plip_nl_pistack,plip_nl_cation,plip_nl_salt,plip_nl_metal,plip_halogen,plip_water,plip_nl_bb,plip_nl_total," +
                                "energy,egor," +
                                "asa,loopFraction," +
                                "plip_betweenness,plip_closeness,plip_clusteringcoefficient," +
                                "plip_hydrogen_betweenness,plip_hydrogen_closeness,plip_hydrogen_clusteringcoefficient," +
                                "plip_hydrophobic_betweenness,plip_hydrophobic_closeness,plip_hydrophobic_clusteringcoefficient," +
                                "conv_total,conv_l_total,conv_nl_total," +
                                "conv_betweenness,conv_closeness,conv_clusteringcoefficient," +
                                "plip_distinct_neighborhoods,conv_distinct_neighborhoods," +
                                "misfolded" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.DATA_DIRECTORY.resolve("native")
                        .resolve("statistics")
                        .resolve("misfolded.csv"),
                output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");

            return Optional.of(handleProtein(split[0].split("_")[0], "misfolded").get() + System.lineSeparator() +
                    handleProtein(split[1].split("_")[0], "native").get());
        } catch (Exception e) {
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }

    private static Optional<String> handleProtein(String pdbId, String misfolded) {
        System.out.println("handling " + misfolded + " protein " + pdbId);
        Structure structure = StructureParser.fromPdbId(pdbId).parse();
        return handleChain(pdbId, structure.getFirstChain(), misfolded);
    }

    private static Optional<String> handleChain(String pdbId, Chain chain, String misfolded) {
        ResidueGraph conventionalProteinGraph = ResidueGraph.createDistanceResidueGraph(chain);

        if(misfolded.equals("misfolded")) {
            Path xmlPath = Paths.get("/home/bittrich/git/phd_sb_repo/data/native/xml/" + pdbId + ".xml");
            if(!Files.exists(xmlPath)) {
                return Optional.empty();
            }
            try {
                PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                        Jsoup.parse(xmlPath.toFile(),
                                "UTF-8"));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        return Optional.of(chain.aminoAcids()
                .map(aminoAcid -> {
                    GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

                    PLIPInteractionContainer plipInteractionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                    PLIPInteractionContainer nonLocalPlipInteractionContainer = new PLIPInteractionContainer(null,
                            plipInteractionContainer
                                    .getInteractions()
                                    .stream()
                                    // interactions have to be non-local
                                    .filter(inter -> Math.abs(inter.getPartner1().getResidueIndex() - inter.getPartner2().getResidueIndex()) > 5)
                                    .collect(Collectors.toList()));
                    PLIPInteractionContainer localPlipInteractionContainer = new PLIPInteractionContainer(null,
                            plipInteractionContainer
                                    .getInteractions()
                                    .stream()
                                    // interactions have to be local
                                    .filter(inter -> !nonLocalPlipInteractionContainer.getInteractions().contains(inter))
                                    .collect(Collectors.toList()));

                    ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer = aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);

                    return pdbId + "," +
                            "A" + "," +
                            aminoAcid.getResidueIdentifier() + "," +
                            aminoAcid.getOneLetterCode() + "," +
                            sse.getSecondaryStructure().getReducedRepresentation() + "," +
                            sse.getSecondaryStructure().getOneLetterRepresentation() + "," +
                            sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize() + "," +
                            (aminoAcid.getFeature(AccessibleSurfaceArea.class).isExposed() ? "exposed" : "buried") + "," +

                            plipInteractionContainer.getHydrogenBonds().size() + "," +
                            plipInteractionContainer.getHydrophobicInteractions().size() + "," +
                            plipInteractionContainer.getPiStackings().size() + "," +
                            plipInteractionContainer.getPiCationInteractions().size() + "," +
                            plipInteractionContainer.getSaltBridges().size() + "," +
                            plipInteractionContainer.getMetalComplexes().size() + "," +
                            plipInteractionContainer.getHalogenBonds().size() + "," +
                            plipInteractionContainer.getWaterBridges().size() + "," +
                            plipInteractionContainer.getBackboneInteractions().size() + "," +
                            plipInteractionContainer.getInteractions().size() + "," +

                            localPlipInteractionContainer.getHydrogenBonds().size() + "," +
                            localPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
                            localPlipInteractionContainer.getPiStackings().size() + "," +
                            localPlipInteractionContainer.getPiCationInteractions().size() + "," +
                            localPlipInteractionContainer.getSaltBridges().size() + "," +
                            localPlipInteractionContainer.getMetalComplexes().size() + "," +
                            localPlipInteractionContainer.getHalogenBonds().size() + "," +
                            localPlipInteractionContainer.getWaterBridges().size() + "," +
                            localPlipInteractionContainer.getBackboneInteractions().size() + "," +
                            localPlipInteractionContainer.getInteractions().size() + "," +

                            nonLocalPlipInteractionContainer.getHydrogenBonds().size() + "," +
                            nonLocalPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
                            nonLocalPlipInteractionContainer.getPiStackings().size() + "," +
                            nonLocalPlipInteractionContainer.getPiCationInteractions().size() + "," +
                            nonLocalPlipInteractionContainer.getSaltBridges().size() + "," +
                            nonLocalPlipInteractionContainer.getMetalComplexes().size() + "," +
                            nonLocalPlipInteractionContainer.getHalogenBonds().size() + "," +
                            nonLocalPlipInteractionContainer.getWaterBridges().size() + "," +
                            nonLocalPlipInteractionContainer.getBackboneInteractions().size() + "," +
                            nonLocalPlipInteractionContainer.getInteractions().size() + "," +

                            StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                            StandardFormat.format(aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction()) + "," +

                            StandardFormat.format(aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                            StandardFormat.format(aminoAcid.getFeature(LoopFraction.class).getLoopFraction()) + "," +

                            StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getBetweenness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getCloseness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getClusteringCoefficient()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getHydrogenPlip().getBetweenness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getHydrogenPlip().getCloseness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getHydrogenPlip().getClusteringCoefficient()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getHydrophobicPlip().getBetweenness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getHydrophobicPlip().getCloseness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getHydrophobicPlip().getClusteringCoefficient()) + "," +

                            conventionalProteinGraph.getContactsOf(aminoAcid).size() + "," +
                            conventionalProteinGraph.getLocalContactsOf(aminoAcid).size() + "," +
                            conventionalProteinGraph.getNonLocalContactsOf(aminoAcid).size() + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getBetweenness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getCloseness()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getClusteringCoefficient()) + "," +

                            StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getDistinctNeighborhoodCount()) + "," +
                            StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getDistinctNeighborhoodCount()) + "," +

                            misfolded;
                })
                .collect(Collectors.joining(System.lineSeparator())));
    }
}

package de.bioforscher.thermostability.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.graphs.ProteinGraph;
import de.bioforscher.jstructure.feature.graphs.ProteinGraphFactory;
import de.bioforscher.jstructure.feature.graphs.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.start2fold.Start2FoldConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A02_WriteDatasetCsv {
    private static final Logger logger = LoggerFactory.getLogger(A02_WriteDatasetCsv.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.THERMOSTABILITY_DIRECTORY.resolve("sadeghi.list"))
                .skip(37)
                .limit(20)
                .map(A02_WriteDatasetCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "exposed," +
                                "plip_hydrogen,plip_hydrophobic,plip_bb,plip_total," +
                                "plip_l_hydrogen,plip_l_hydrophobic,plip_l_bb,plip_l_total," +
                                "plip_nl_hydrogen,plip_nl_hydrophobic,plip_nl_bb,plip_nl_total," +
                                "energy,egor," +
                                "asa,loopFraction," +
                                "plip_betweenness,plip_closeness,plip_clusteringcoefficient," +
                                "plip_hydrogen_betweenness,plip_hydrogen_closeness,plip_hydrogen_clusteringcoefficient," +
                                "plip_hydrophobic_betweenness,plip_hydrophobic_closeness,plip_hydrophobic_clusteringcoefficient," +
                                "conv_total,conv_l_total,conv_nl_total," +
                                "conv_betweenness,conv_closeness,conv_clusteringcoefficient," +
                                "plip_distinct_neighborhoods,conv_distinct_neighborhoods," +
                                "thermophile" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.THERMOSTABILITY_DIRECTORY.resolve("statistics").resolve("thermophile-2.csv"),
                output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");

            return Optional.of(handleProtein(split[0], split[1], "thermophile") + System.lineSeparator() +
                    handleProtein(split[2], split[3], "mesophile"));
        } catch (Exception e) {
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }

    private static String handleProtein(String pdbId, String chainIds, String thermophile) {
        System.out.println("handling " + thermophile + " protein " + pdbId + " with chains " + chainIds);
        Structure structure = StructureParser.fromPdbId(pdbId).parse();
        List<Chain> chains = chainIds.equals("Null") ?
                Stream.of(structure.chainsWithAminoAcids().findFirst().get()).collect(Collectors.toList()) :
                Pattern.compile(",").splitAsStream(chainIds)
                        .map(chainId -> {
                            Optional<Chain> optionalChain = structure.select().chainId(chainId).asOptionalChain();
                            return optionalChain.orElseGet(() -> structure.select().chainId(chainId.toUpperCase()).asChain());
                        })
                        .collect(Collectors.toList());
        return chains.stream()
                .map(chain -> handleChain(pdbId, chain, thermophile))
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static String handleChain(String pdbId, Chain chain, String thermophile) {
        ProteinGraph conventionalProteinGraph = ProteinGraphFactory.createProteinGraph(chain, ProteinGraphFactory.InteractionScheme.CALPHA8);

        return chain.aminoAcids()
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
                            plipInteractionContainer.getBackboneInteractions().size() + "," +
                            plipInteractionContainer.getInteractions().size() + "," +

                            localPlipInteractionContainer.getHydrogenBonds().size() + "," +
                            localPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
                            localPlipInteractionContainer.getBackboneInteractions().size() + "," +
                            localPlipInteractionContainer.getInteractions().size() + "," +

                            nonLocalPlipInteractionContainer.getHydrogenBonds().size() + "," +
                            nonLocalPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
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

                            thermophile;
                })
                .collect(Collectors.joining(System.lineSeparator()));
    }
}

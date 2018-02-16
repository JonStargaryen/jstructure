package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.graphs.*;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.interactions.PLIPRestServiceQuery;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A04_WriteTransitionStateCsv {
    //TODO add plip backbone contacts
    private static final Logger logger = LoggerFactory.getLogger(A04_WriteTransitionStateCsv.class);
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) throws IOException {
        String localOutput = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A04_WriteTransitionStateCsv::handleLineLocally)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa," +

                                "plip_total," +
                                "plip_l_total," +
                                "plip_nl_total," +
                                "plip_betweenness,plip_closeness,plip_clusteringcoefficient," +
                                "plip_hydrogen_total," +
                                "plip_hydrogen_l_total," +
                                "plip_hydrogen_nl_total," +
                                "plip_hydrogen_betweenness,plip_hydrogen_closeness,plip_hydrogen_clusteringcoefficient," +
                                "plip_hydrophobic_total," +
                                "plip_hydrophobic_l_total," +
                                "plip_hydrophobic_nl_total," +
                                "plip_hydrophobic_betweenness,plip_hydrophobic_closeness,plip_hydrophobic_clusteringcoefficient," +
                                "conv_total," +
                                "conv_l_total," +
                                "conv_nl_total," +
                                "conv_betweenness,conv_closeness,conv_clusteringcoefficient," +
                                "plip_distinct_neighborhoods," +
                                "conv_distinct_neighborhoods," +
                                "energy," +
                                "asa,loopFraction," +
                                "folds,transition" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("transition-state.csv"),
                localOutput);
    }

    private static Optional<String> handleLineLocally(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

            Structure structure = StructureParser.source(pdbId).parse();
            Chain originalChain = structure.chains().findFirst().get();
            ProteinGraph originalFullPlipGraph = ProteinGraphFactory.createProteinGraph(originalChain,
                    ProteinGraphFactory.InteractionScheme.SALENTIN2015);
            ProteinGraph originalHydrogenPlipGraph = ProteinGraphFactory.createProteinGraph(originalChain,
                    ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROGEN_BONDS);
            ProteinGraph originalHydrophobicPlipGraph = ProteinGraphFactory.createProteinGraph(originalChain,
                    ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROPHOBIC_INTERACTION);
            ProteinGraph originalConvGraph = ProteinGraphFactory.createProteinGraph(originalChain,
                    ProteinGraphFactory.InteractionScheme.CALPHA8);

            Start2FoldXmlParser.parseSpecificExperiment(originalChain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            List<AminoAcid> earlyFoldingResidues = originalChain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            List<Chain> reconstructedChains = Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/" +
                    "reconstruction-start2fold/reconstructions/" + pdbId + "-early-conventional-1/stage1/"))
                    .filter(path -> path.toFile().getName().contains("_model"))
                    .map(path -> StructureParser.source(path)
                            .forceProteinName(IdentifierFactory.createProteinIdentifier(pdbId, path.toFile().getName().split("_")[2].split("\\.")[0]))
                            .parse()
                            .getChains()
                            .get(0))
                    .collect(Collectors.toList());

            for (Chain reconstructedChain : reconstructedChains) {
                Document document = PLIPRestServiceQuery.calculateIntraChainDocument(reconstructedChain);
                PLIP_INTRA_MOLECULAR_ANNOTATOR.process(originalChain, document);
            }

            List<ProteinGraph> convGraphs = reconstructedChains.stream()
                    .map(c -> ProteinGraphFactory.createProteinGraph(c, ProteinGraphFactory.InteractionScheme.CALPHA8))
                    .collect(Collectors.toList());
            List<ProteinGraphCalculations> convGraphCalculations = convGraphs.stream()
                    .map(ProteinGraphCalculations::new)
                    .collect(Collectors.toList());
            List<ProteinGraph> fullPlipGraphs = reconstructedChains.stream()
                    .map(c -> ProteinGraphFactory.createProteinGraph(c, ProteinGraphFactory.InteractionScheme.SALENTIN2015))
                    .collect(Collectors.toList());
            List<ProteinGraphCalculations> fullPlipGraphCalculations = fullPlipGraphs.stream()
                    .map(ProteinGraphCalculations::new)
                    .collect(Collectors.toList());
            List<ProteinGraph> hydrogenPlipGraphs = reconstructedChains.stream()
                    .map(c -> ProteinGraphFactory.createProteinGraph(c, ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROGEN_BONDS))
                    .collect(Collectors.toList());
            List<ProteinGraphCalculations> hydrogenPlipGraphCalculations = fullPlipGraphs.stream()
                    .map(ProteinGraphCalculations::new)
                    .collect(Collectors.toList());
            List<ProteinGraph> hydrophobicPlipGraphs = reconstructedChains.stream()
                    .map(c -> ProteinGraphFactory.createProteinGraph(c, ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROPHOBIC_INTERACTION))
                    .collect(Collectors.toList());
            List<ProteinGraphCalculations> hydrophobicPlipGraphCalculations = fullPlipGraphs.stream()
                    .map(ProteinGraphCalculations::new)
                    .collect(Collectors.toList());

            return Optional.of(originalChain.aminoAcids()
                    .map(aminoAcid -> {
                        ResidueTopologicPropertiesContainer container = aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);
                        ResidueIdentifier residueIdentifier = aminoAcid.getResidueIdentifier();

                        return pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +

                                originalFullPlipGraph.getContactsOf(aminoAcid).size() + "," +
                                originalFullPlipGraph.getLocalContactsOf(aminoAcid).size() + "," +
                                originalFullPlipGraph.getNonLocalContactsOf(aminoAcid).size() + "," +
                                StandardFormat.format(container.getFullPlip().getBetweenness()) + "," +
                                StandardFormat.format(container.getFullPlip().getCloseness()) + "," +
                                StandardFormat.format(container.getFullPlip().getClusteringCoefficient()) + "," +

                                originalHydrogenPlipGraph.getContactsOf(aminoAcid).size() + "," +
                                originalHydrogenPlipGraph.getLocalContactsOf(aminoAcid).size() + "," +
                                originalHydrogenPlipGraph.getNonLocalContactsOf(aminoAcid).size() + "," +
                                StandardFormat.format(container.getHydrogenPlip().getBetweenness()) + "," +
                                StandardFormat.format(container.getHydrogenPlip().getCloseness()) + "," +
                                StandardFormat.format(container.getHydrogenPlip().getClusteringCoefficient()) + "," +

                                originalHydrophobicPlipGraph.getContactsOf(aminoAcid).size() + "," +
                                originalHydrophobicPlipGraph.getLocalContactsOf(aminoAcid).size() + "," +
                                originalHydrophobicPlipGraph.getNonLocalContactsOf(aminoAcid).size() + "," +
                                StandardFormat.format(container.getHydrophobicPlip().getBetweenness()) + "," +
                                StandardFormat.format(container.getHydrophobicPlip().getCloseness()) + "," +
                                StandardFormat.format(container.getHydrophobicPlip().getClusteringCoefficient()) + "," +

                                originalConvGraph.getContactsOf(aminoAcid).size() + "," +
                                originalConvGraph.getLocalContactsOf(aminoAcid).size() + "," +
                                originalConvGraph.getNonLocalContactsOf(aminoAcid).size() + "," +
                                StandardFormat.format(container.getConventional().getBetweenness()) + "," +
                                StandardFormat.format(container.getConventional().getCloseness()) + "," +
                                StandardFormat.format(container.getConventional().getClusteringCoefficient()) + "," +

                                container.getFullPlip().getDistinctNeighborhoodCount() + "," +
                                container.getConventional().getDistinctNeighborhoodCount() + "," +

                                StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +

                                StandardFormat.format(aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(LoopFraction.class).getLoopFraction()) + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                                "native" + System.lineSeparator() +

                                pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +

                                StandardFormat.format(fullPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(fullPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(fullPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getNonLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(fullPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.betweenness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(fullPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.closeness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(fullPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.clusteringCoefficient(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +

                                StandardFormat.format(hydrogenPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrogenPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrogenPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getNonLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrogenPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.betweenness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrogenPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.closeness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrogenPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.clusteringCoefficient(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +

                                StandardFormat.format(hydrophobicPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrophobicPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrophobicPlipGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getNonLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrophobicPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.betweenness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrophobicPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.closeness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(hydrophobicPlipGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.clusteringCoefficient(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +

                                StandardFormat.format(convGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(convGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(convGraphs.stream()
                                        .mapToInt(proteinGraph -> proteinGraph.getNonLocalContactsOf(residueIdentifier).size())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(convGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.betweenness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(convGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.closeness(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(convGraphCalculations.stream()
                                        .mapToDouble(proteinGraphCalculations -> proteinGraphCalculations.clusteringCoefficient(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +

                                StandardFormat.format(fullPlipGraphCalculations.stream()
                                        .mapToInt(proteinGraphCalculations -> proteinGraphCalculations.distinctNeighborhoodCount(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(convGraphCalculations.stream()
                                        .mapToInt(proteinGraphCalculations -> proteinGraphCalculations.distinctNeighborhoodCount(residueIdentifier))
                                        .average()
                                        .getAsDouble()) + "," +

                                StandardFormat.format(reconstructedChains.stream()
                                        .map(chain -> chain.select().residueIdentifier(aminoAcid.getResidueIdentifier()).asAminoAcid())
                                        .mapToDouble(aa -> aa.getFeature(EnergyProfile.class).getSolvationEnergy())
                                        .average()
                                        .getAsDouble()) + "," +

                                StandardFormat.format(reconstructedChains.stream()
                                        .map(chain -> chain.select().residueIdentifier(aminoAcid.getResidueIdentifier()).asAminoAcid())
                                        .mapToDouble(aa -> aa.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea())
                                        .average()
                                        .getAsDouble()) + "," +
                                StandardFormat.format(reconstructedChains.stream()
                                        .map(chain -> chain.select().residueIdentifier(aminoAcid.getResidueIdentifier()).asAminoAcid())
                                        .mapToDouble(aa -> aa.getFeature(LoopFraction.class).getLoopFraction())
                                        .average()
                                        .getAsDouble()) + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                                "transition";
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }
}

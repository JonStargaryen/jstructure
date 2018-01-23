package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.graphs.ProteinGraphFactory;
import de.bioforscher.jstructure.feature.graphs.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A04_WriteTransitionStateCsv {
    private static final Logger logger = LoggerFactory.getLogger(A04_WriteTransitionStateCsv.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A04_WriteTransitionStateCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa," +
                                "conv_total," +
                                "conv_l_total," +
                                "conv_nl_total," +
                                "conv_closeness,conv_clusteringcoefficient," +
                                "conv_distinct_neighborhoods," +
                                "folds,transition" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("transition-state.csv"),
                output);
    }

    private static Optional<String> handleLine(String line) {
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
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            Graph<AminoAcid> proteinGraph = ProteinGraphFactory.createProteinGraph(chain,
                    ProteinGraphFactory.InteractionScheme.CALPHA8);
            List<Edge<AminoAcid>> interactions = proteinGraph.getEdges();
            List<Edge<AminoAcid>> localInteractions = interactions.stream()
                    .filter(edge -> Math.abs(edge.getLeft().getResidueIndex() - edge.getRight().getResidueIndex()) <= 6)
                    .collect(Collectors.toList());
            List<Edge<AminoAcid>> nonLocalInteractions = interactions.stream()
                    .filter(edge -> Math.abs(edge.getLeft().getResidueIndex() - edge.getRight().getResidueIndex()) > 6)
                    .collect(Collectors.toList());

            List<Chain> reconstructedChains = Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/" +
                    "reconstruction-start2fold/reconstructions/" + pdbId + "-early-conventional-1/stage1/"))
                    .filter(path -> path.toFile().getName().contains("_model"))
                    .map(path -> StructureParser.source(path).parse().getChains().get(0))
                    .collect(Collectors.toList());
            List<Graph<AminoAcid>> reconstructedGraphs = reconstructedChains.stream()
                    .map(c -> ProteinGraphFactory.createProteinGraph(c, ProteinGraphFactory.InteractionScheme.CALPHA8))
                    .collect(Collectors.toList());

            return Optional.of(chain.aminoAcids()
                    .map(aminoAcid -> {
                        ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer = aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);
                        int distinctNeighborhoodCount = determineDistinctNeighborhoodCount(interactions, aminoAcid);

                        double averageInteractionCount = reconstructedGraphs.stream()
                                .mapToInt(graph -> filterInteractions(graph, aminoAcid))
                                .average()
                                .getAsDouble();
                        double averageLocalInteractionCount = reconstructedGraphs.stream()
                                .mapToInt(graph -> filterLocalInteractions(graph, aminoAcid))
                                .average()
                                .getAsDouble();
                        double averageNonLocalInteractionCount = reconstructedGraphs.stream()
                                .mapToInt(graph -> filterNonLocalInteractions(graph, aminoAcid))
                                .average()
                                .getAsDouble();
                        double averageCloseness = reconstructedGraphs.stream()
                                .mapToDouble(graph -> {
                                    AminoAcid aa = graph.nodes().filter(node -> node.getResidueIdentifier().equals(aminoAcid.getResidueIdentifier())).findFirst().get();
                                    return graph.calculate().closeness(aa);
                                })
                                .average()
                                .getAsDouble();
                        double averageClusteringCoefficient = reconstructedGraphs.stream()
                                .mapToDouble(graph -> {
                                    AminoAcid aa = graph.nodes().filter(node -> node.getResidueIdentifier().equals(aminoAcid.getResidueIdentifier())).findFirst().get();
                                    return graph.calculate().clusteringCoefficient(aa);
                                })
                                .average()
                                .getAsDouble();
                        double averageDistinctNeighborhoodCount = reconstructedGraphs.stream()
                                .mapToInt(graph -> {
                                    AminoAcid aa = graph.nodes().filter(node -> node.getResidueIdentifier().equals(aminoAcid.getResidueIdentifier())).findFirst().get();
                                    return determineDistinctNeighborhoodCount(graph.getEdges(), aa);
                                })
                                .average()
                                .getAsDouble();

                        return pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +

                                interactions.stream()
                                        .filter(edge -> edge.contains(aminoAcid)).count() + "," +
                                localInteractions.stream()
                                        .filter(edge -> edge.contains(aminoAcid)).count() + "," +
                                nonLocalInteractions.stream()
                                        .filter(edge -> edge.contains(aminoAcid)).count() + "," +

                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getCloseness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getClusteringCoefficient()) + "," +

                                distinctNeighborhoodCount + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                                "native" + System.lineSeparator() +

                                pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +

                                StandardFormat.format(averageInteractionCount) + "," +
                                StandardFormat.format(averageLocalInteractionCount) + "," +
                                StandardFormat.format(averageNonLocalInteractionCount) + "," +

                                StandardFormat.format(averageCloseness) + "," +
                                StandardFormat.format(averageClusteringCoefficient) + "," +

                                StandardFormat.format(averageDistinctNeighborhoodCount) + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                                "transition";
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }

    private static int determineDistinctNeighborhoodCount(List<Edge<AminoAcid>> interactions, AminoAcid aminoAcid) {
        List<Integer> interactingGroups = interactions.stream()
                .filter(plipInteraction -> aminoAcid.equals(plipInteraction.getLeft()) || aminoAcid.equals(plipInteraction.getRight()))
                .filter(inter -> Math.abs(inter.getLeft().getResidueIndex() - inter.getRight().getResidueIndex()) > 6)
                .map(plipInteraction -> aminoAcid.equals(plipInteraction.getLeft()) ? plipInteraction.getRight() : plipInteraction.getLeft())
                .map(Group::getResidueIndex)
                .collect(Collectors.toList());
        List<Integer> distinctNeighborhoods = new ArrayList<>();
        for(int interactingGroup : interactingGroups) {
            if(distinctNeighborhoods.stream()
                    .noneMatch(residueIndex -> Math.abs(residueIndex - interactingGroup) > 6)) {
                distinctNeighborhoods.add(interactingGroup);
            }
        }
        return distinctNeighborhoods.size();
    }

    private static int filterInteractions(Graph<AminoAcid> graph, AminoAcid aminoAcid) {
        return (int) graph.edges()
                .filter(edge -> filterEdge(edge, aminoAcid))
                .count();
    }

    private static int filterLocalInteractions(Graph<AminoAcid> graph, AminoAcid aminoAcid) {
        return (int) graph.edges()
                .filter(edge -> Math.abs(edge.getLeft().getResidueIndex() - edge.getRight().getResidueIndex()) <= 6)
                .filter(edge -> filterEdge(edge, aminoAcid))
                .count();
    }

    private static int filterNonLocalInteractions(Graph<AminoAcid> graph, AminoAcid aminoAcid) {
        return (int) graph.edges()
                .filter(edge -> Math.abs(edge.getLeft().getResidueIndex() - edge.getRight().getResidueIndex()) > 6)
                .filter(edge -> filterEdge(edge, aminoAcid))
                .count();
    }

    private static boolean filterEdge(Edge<AminoAcid> edge, AminoAcid aminoAcid) {
        return edge.getLeft().getResidueIdentifier().equals(aminoAcid.getResidueIdentifier()) ||
                edge.getRight().getResidueIdentifier().equals(aminoAcid.getResidueIdentifier());
    }
}

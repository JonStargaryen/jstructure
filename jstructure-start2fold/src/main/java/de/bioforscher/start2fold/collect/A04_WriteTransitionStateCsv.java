package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.feature.graphs.ProteinGraph;
import de.bioforscher.jstructure.feature.graphs.ProteinGraphFactory;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Chain;
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

            ProteinGraph proteinGraph = ProteinGraphFactory.createProteinGraph(chain,
                    ProteinGraphFactory.InteractionScheme.CALPHA8);
            List<Pair<AminoAcid, AminoAcid>> contacts = proteinGraph.getContacts();
            List<Pair<AminoAcid, AminoAcid>> localInteractions = proteinGraph.getLocalContacts();
            List<Pair<AminoAcid, AminoAcid>> nonLocalInteractions = proteinGraph.getNonLocalContacts();

            List<Chain> reconstructedChains = Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/" +
                    "reconstruction-start2fold/reconstructions/" + pdbId + "-early-conventional-1/stage1/"))
                    .filter(path -> path.toFile().getName().contains("_model"))
                    .map(path -> StructureParser.source(path).parse().getChains().get(0))
                    .collect(Collectors.toList());

            //TODO compute PLIP
            List<ProteinGraph> reconstructedGraphs = reconstructedChains.stream()
                    .map(c -> ProteinGraphFactory.createProteinGraph(c, ProteinGraphFactory.InteractionScheme.CALPHA8))
                    .collect(Collectors.toList());

            //TODO reimpl
//            return Optional.of(chain.aminoAcids()
//                    .map(aminoAcid -> {
//                        ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer = aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);
//                        int distinctNeighborhoodCount = determineDistinctNeighborhoodCount(interactions, aminoAcid);
//
//                        double averageInteractionCount = reconstructedGraphs.stream()
//                                .mapToInt(graph -> filterInteractions(graph, aminoAcid))
//                                .average()
//                                .getAsDouble();
//                        double averageLocalInteractionCount = reconstructedGraphs.stream()
//                                .mapToInt(graph -> filterLocalInteractions(graph, aminoAcid))
//                                .average()
//                                .getAsDouble();
//                        double averageNonLocalInteractionCount = reconstructedGraphs.stream()
//                                .mapToInt(graph -> filterNonLocalInteractions(graph, aminoAcid))
//                                .average()
//                                .getAsDouble();
//                        double averageCloseness = reconstructedGraphs.stream()
//                                .mapToDouble(graph -> {
//                                    AminoAcid aa = graph.nodes().filter(node -> node.getResidueIdentifier().equals(aminoAcid.getResidueIdentifier())).findFirst().get();
//                                    return graph.calculate().closeness(aa);
//                                })
//                                .average()
//                                .getAsDouble();
//                        double averageClusteringCoefficient = reconstructedGraphs.stream()
//                                .mapToDouble(graph -> {
//                                    AminoAcid aa = graph.nodes().filter(node -> node.getResidueIdentifier().equals(aminoAcid.getResidueIdentifier())).findFirst().get();
//                                    return graph.calculate().clusteringCoefficient(aa);
//                                })
//                                .average()
//                                .getAsDouble();
//                        double averageDistinctNeighborhoodCount = reconstructedGraphs.stream()
//                                .mapToInt(graph -> {
//                                    AminoAcid aa = graph.nodes().filter(node -> node.getResidueIdentifier().equals(aminoAcid.getResidueIdentifier())).findFirst().get();
//                                    return determineDistinctNeighborhoodCount(graph.getEdges(), aa);
//                                })
//                                .average()
//                                .getAsDouble();
//
//                        return pdbId + "," +
//                                "A" + "," +
//                                aminoAcid.getResidueIdentifier() + "," +
//                                aminoAcid.getOneLetterCode() + "," +
//
//                                contacts.stream()
//                                        .filter(edge -> edge.contains(aminoAcid)).count() + "," +
//                                localInteractions.stream()
//                                        .filter(edge -> edge.contains(aminoAcid)).count() + "," +
//                                nonLocalInteractions.stream()
//                                        .filter(edge -> edge.contains(aminoAcid)).count() + "," +
//
//                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getCloseness()) + "," +
//                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getClusteringCoefficient()) + "," +
//
//                                distinctNeighborhoodCount + "," +
//
//                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
//                                "native" + System.lineSeparator() +
//
//                                pdbId + "," +
//                                "A" + "," +
//                                aminoAcid.getResidueIdentifier() + "," +
//                                aminoAcid.getOneLetterCode() + "," +
//
//                                StandardFormat.format(averageInteractionCount) + "," +
//                                StandardFormat.format(averageLocalInteractionCount) + "," +
//                                StandardFormat.format(averageNonLocalInteractionCount) + "," +
//
//                                StandardFormat.format(averageCloseness) + "," +
//                                StandardFormat.format(averageClusteringCoefficient) + "," +
//
//                                StandardFormat.format(averageDistinctNeighborhoodCount) + "," +
//
//                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
//                                "transition";
//                    })
//                    .collect(Collectors.joining(System.lineSeparator())));
            return Optional.empty();
        } catch (Exception e) {
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }
}

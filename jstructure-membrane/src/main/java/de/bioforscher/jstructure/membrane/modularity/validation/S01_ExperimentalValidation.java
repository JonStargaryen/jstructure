package de.bioforscher.jstructure.membrane.modularity.validation;

import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.clustering.ClusteringScorer;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.modularity.GraphFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class S01_ExperimentalValidation {
    private static final Logger logger = LoggerFactory.getLogger(S01_ExperimentalValidation.class);

    public static void main(String[] args) {
        Path dataset = MembraneConstants.MODULARITY_DATASET_DIRECTORY;

        String output = MembraneConstants.lines(dataset.resolve("ids.list"))
                .limit(1)
                .map(id -> handleId(dataset, id))
                .collect(Collectors.joining(System.lineSeparator(),
                        "size,mod_num_exp," +
                                "mod_num_sse,sse_fm," +
                                getInflationValueRange()
                                        .mapToObj(inflation -> "wmcl_" + inflation)
                                        .map(inflation -> "mod_num_" + inflation + "," + inflation + "_fm")
                                        .collect(Collectors.joining(",")) + "," +
                                getInflationValueRange()
                                        .mapToObj(inflation -> "uwmcl_" + inflation)
                                        .map(inflation -> "mod_num_" + inflation + "," + inflation + "_fm")
                                        .collect(Collectors.joining(",")) + System.lineSeparator(),
                        ""));

        System.out.println(output);
    }

    private static String handleId(Path dataset, String id) {
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];

        logger.info("handling chain {} of {}",
                chainId,
                pdbId);

        Chain chain = StructureParser.source(dataset.resolve("pdb").resolve(pdbId + ".pdb"))
                .minimalParsing(true)
                .parse()
                .select()
                .chainName(chainId)
                .asChain();

        PartitionedGraph<AminoAcid> experimental = GraphFactory.createPartitionedGraphFromDefinitionFile(chain,
                dataset.resolve("exp").resolve(id + ".exp"));
        PartitionedGraph<AminoAcid> secondaryStructure = GraphFactory.createPartitionedGraphFromSecondaryStructureElements(chain);
        System.out.println(experimental.getDefinitionString());
        System.out.println(secondaryStructure.getDefinitionString());
        Map<Integer, PartitionedGraph<AminoAcid>> weightedMclResults = getInflationValueRange()
                .boxed()
                .collect(Collectors.toMap(Function.identity(),
                        inflation -> {
                            double i = inflation / (double) 10;
                            logger.info("sampling graph of {} at inflation {}",
                                    chain.getChainIdentifier(),
                                    i);
                            return GraphFactory.createPartionedGraphFromWeightedPlipData(chain,
                                    dataset.resolve("plip").resolve(id + ".plip"),
                                    i);
                        }));
        Map<Integer, PartitionedGraph<AminoAcid>> unweightedMclResults = getInflationValueRange()
                .boxed()
                .collect(Collectors.toMap(Function.identity(),
                        inflation -> {
                            double i = inflation / (double) 10;
                            logger.info("sampling graph of {} at inflation {}",
                                    chain.getChainIdentifier(),
                                    i);
                            return GraphFactory.createPartitionedGraphFromUnweightedPlipData(chain,
                                    dataset.resolve("plip").resolve(id + ".plip"),
                                    i);
                        }));

        return chain.aminoAcids().count() + "," + experimental.getNumberOfModules() + "," +
                composeScoreString(experimental, secondaryStructure) + "," +
                getInflationValueRange()
                        .mapToObj(weightedMclResults::get)
                        .map(mclResult -> composeScoreString(experimental, mclResult))
                        .collect(Collectors.joining(",")) + "," +
                getInflationValueRange()
                        .mapToObj(unweightedMclResults::get)
                        .map(mclResult -> composeScoreString(experimental, mclResult))
                        .collect(Collectors.joining(","));
    }

    private static IntStream getInflationValueRange() {
        // originally, values are in the interval [1.2,5.0] - big values lead to fragmented clusters - default is 2.0
        return IntStream.range(12, 21);
    }

    private static String composeScoreString(PartitionedGraph<AminoAcid> reference, PartitionedGraph<AminoAcid> clustering) {
        return clustering.getNumberOfModules() + "," +
                ClusteringScorer.fowlkesMallowsIndex(reference, clustering) /*+ "," +
                ClusteringScorer.naiveScore(reference, clustering) + "," +
                Math.log(ClusteringScorer.chiSquaredCoefficient(reference, clustering))*/;
    }
}

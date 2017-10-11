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

/**
 * Compare experimental, NetCarto and MCL results by varying scores.
 */
public class CommunityDetectionComparison {
    private static final Logger logger = LoggerFactory.getLogger(CommunityDetectionComparison.class);

    public static void main(String[] args) {
        Path dataset = MembraneConstants.MODULARITY_DATASET_DIRECTORY;

        String output = MembraneConstants.lines(dataset.resolve("ids.list"))
                .map(id -> handleId(dataset, id))
                .collect(Collectors.joining(System.lineSeparator(),
                        "num_exp," +
                                "num_netcarto,netcarto_randindex,netcarto_naivescore,netcarto_logchisquared," +
                                getInflationValueRange()
                                        .mapToObj(inflation -> "mcl_" + inflation)
                                        .map(inflation -> "num_" + inflation + "," + inflation + "_randindex," + inflation + "_naivescore," + inflation + "_logchisquared")
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
        PartitionedGraph<AminoAcid> netCarto = GraphFactory.createPartitionedGraphFromNetCartoFile(chain,
                dataset.resolve("network").resolve(id + "_plip.modules.dat"));
        Map<Integer, PartitionedGraph<AminoAcid>> mclResults = getInflationValueRange()
                .boxed()
                .collect(Collectors.toMap(Function.identity(),
                        inflation -> {
                            double i = inflation / (double) 10;
                            logger.info("sampling graph of {} at inflation {}",
                                    chain.getChainIdentifier(),
                                    i);
                            return GraphFactory.createPartitionedGraphFromMCL(chain,
                                    dataset.resolve("plip").resolve(id + ".plip"),
                                    i);
                        }));

        return experimental.getNumberOfModules() + "," +
                composeScoreString(experimental, netCarto) + "," +
                getInflationValueRange()
                        .mapToObj(mclResults::get)
                        .map(mclResult -> composeScoreString(experimental, mclResult))
                        .collect(Collectors.joining(","));
    }

    private static IntStream getInflationValueRange() {
        // originally, values are in the interval [1.2,5.0] - big values lead to fragmented clusters - default is 2.0
        return IntStream.range(12, 21);
    }

    private static String composeScoreString(PartitionedGraph<AminoAcid> partitionedGraph1, PartitionedGraph<AminoAcid> partitionedGraph2) {
        return partitionedGraph2.getNumberOfModules() + "," +
                ClusteringScorer.randIndex(partitionedGraph1, partitionedGraph2) + "," +
                ClusteringScorer.naiveScore(partitionedGraph1, partitionedGraph2) + "," +
                Math.log(ClusteringScorer.chiSquaredCoefficient(partitionedGraph1, partitionedGraph2));
    }
}

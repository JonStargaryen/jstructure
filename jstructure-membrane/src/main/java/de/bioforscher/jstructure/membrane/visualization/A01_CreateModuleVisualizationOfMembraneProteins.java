package de.bioforscher.jstructure.membrane.visualization;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.bioforscher.jstructure.contacts.ProteinGraphFactory;
import de.bioforscher.jstructure.contacts.visualization.GraphVisualizer;
import de.bioforscher.jstructure.contacts.visualization.JsonChain;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.impl.MarkovClusteringAlgorithm;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Deprecated
// change alpha and beta directory when needed - russian!
public class A01_CreateModuleVisualizationOfMembraneProteins {
    private static final List<String> ids = Stream.of("1h6s_1", "3o0r_B")
            .collect(Collectors.toList());

    public static void main(String[] args) {
        ids.forEach(A01_CreateModuleVisualizationOfMembraneProteins::handleId);
    }

    private static void handleId(String id) {
        try {
            System.out.println(id);

            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();

            // force outside annotation with PLIP data
            new PLIPIntraMolecularAnnotator().process(structure);

            Graph<AminoAcid> graph = ProteinGraphFactory.createProteinGraph(chain, ProteinGraphFactory.InteractionScheme.SALENTIN2015);

            Map<String, PartitionedGraph<AminoAcid>> inSilicoData = new TreeMap<>();

            PartitionedGraph<AminoAcid> netCartoPartitioning = ProteinGraphFactory.Partitioned.fromNetCartoFile(graph,
                    MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("network/" + id + "_plip.modules.dat"));
            inSilicoData.put("NetCarto", netCartoPartitioning);
            MembraneConstants.write(MembraneConstants.DATA_DIRECTORY.resolve("visualization")
                            .resolve("modules")
                            .resolve("pymol")
                            .resolve(id + "_netCarto.pml"),
                    GraphVisualizer.getPyMolString(netCartoPartitioning));

            for (int inflation = 120; inflation < 205; inflation = inflation + 5) {
                double i = inflation / (double) 100;
                PartitionedGraph<AminoAcid> partitionedGraph = ProteinGraphFactory.Partitioned.fromMCL(graph,
                        MarkovClusteringAlgorithm.DEFAULT_EXPAND_FACTOR,
                        i);

                inSilicoData.put("cla-" + i, partitionedGraph);

                String pyMolString = GraphVisualizer.getPyMolString(partitionedGraph);
                MembraneConstants.write(MembraneConstants.DATA_DIRECTORY.resolve("visualization")
                                .resolve("modules")
                                .resolve("pymol")
                                .resolve(id + "_" + inflation + ".pml"),
                        pyMolString);
            }

            String json = new ObjectMapper().writeValueAsString(new JsonChain(chain, inSilicoData));

            MembraneConstants.write(MembraneConstants.DATA_DIRECTORY.resolve("visualization")
                            .resolve("modules")
                            .resolve("pymol")
                            .resolve(id + ".json"),
                    json);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

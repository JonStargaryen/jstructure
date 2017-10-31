package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.impl.MarkovClusteringAlgorithm;
import de.bioforscher.jstructure.membrane.GraphFactory;
import de.bioforscher.jstructure.membrane.visualization.GraphVisualizer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
import org.jsoup.Jsoup;
import org.junit.Test;

import java.util.stream.Collectors;

public class PlipGraphFactoryTest {
    @Test
    public void shouldPartitionProtein() {
        Chain chain = StructureParser.source("3o0r")
                .minimalParsing(true)
                .parse()
                .select()
                .chainName("B")
                .asChain();

        String document = TestUtils.getResourceAsStream("plip/3o0r_B.plip")
                .collect(Collectors.joining(System.lineSeparator()));
        new PLIPIntraMolecularAnnotator().process(chain, Jsoup.parse(document));
        Graph<AminoAcid> graph = GraphFactory.createProteinGraph(chain, GraphFactory.InteractionScheme.SALENTIN2015);
        PartitionedGraph<AminoAcid> partitionedGraph = GraphFactory.Partitioned.fromMCL(graph, MarkovClusteringAlgorithm.DEFAULT_EXPAND_FACTOR, 2.0);

        System.out.println(GraphVisualizer.getPyMolString(partitionedGraph));
    }
}
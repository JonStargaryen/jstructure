package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.membrane.modularity.visualization.GraphVisualizer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
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
        PartitionedGraph<AminoAcid> graph = GraphFactory.createPartitionedGraphFromPlipData(chain,
                document,
                2.0,
                GraphFactory.WeightingScheme.UNWEIGHTED);

        System.out.println(GraphVisualizer.getPyMolString(graph));
    }
}
package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.clustering.algorithms.MCL;
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
        Graph<AminoAcid> graph = GraphFactory.createGraphFromPlipDocument(chain, document);

        System.out.println(new MCL().clusterGraph(graph).getPyMolString());
    }
}
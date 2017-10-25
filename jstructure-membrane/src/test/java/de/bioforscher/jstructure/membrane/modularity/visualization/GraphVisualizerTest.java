package de.bioforscher.jstructure.membrane.modularity.visualization;

import com.fasterxml.jackson.core.JsonProcessingException;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.membrane.modularity.GraphFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

public class GraphVisualizerTest {
    @Test
    @Ignore
    public void shouldComposeJsonRepresentationFor1f21_A() throws JsonProcessingException {
        String id = "1f21_A";
        composeJsonRepresentation(id);
    }

    @Test
    @Ignore
    public void shouldComposeJsonRepresentationFor3cyt_I() throws JsonProcessingException {
        String id = "3cyt_I";
        composeJsonRepresentation(id);
    }

    private void composeJsonRepresentation(String id) throws JsonProcessingException {
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];

        Structure structure = StructureParser.source(pdbId)
                .minimalParsing(true)
                .parse();
        Chain chain = structure.select()
                .chainName(chainId)
                .asChain();

        List<String> expLines = TestUtils.getResourceAsLines("exp/1f21_A.exp");
        PartitionedGraph<AminoAcid> exp = GraphFactory.createPartitionedGraphFromDefinitionFile(chain, expLines);
        Map<String, PartitionedGraph<AminoAcid>> inSilicoData = new TreeMap<>();

        inSilicoData.put("netcarto", GraphFactory.createPartitionedGraphFromNetCartoFile(chain, TestUtils.getResourceAsLines("netcarto/1f21_A_plip.modules.dat")));
        String plipDocument = TestUtils.getResourceAsStream("plip/1f21_A.plip")
                .collect(Collectors.joining(System.lineSeparator()));

        for(double inflation = 1.2; inflation <= 2.05; inflation = inflation + 0.05) {
            int i = (int) Math.round(inflation * 100);
            System.out.println(i);
            inSilicoData.put("unw-" + i, GraphFactory.createPartitionedGraphFromPlipData(chain,
                    plipDocument,
                    inflation,
                    GraphFactory.WeightingScheme.UNWEIGHTED));
            inSilicoData.put("cla-" + i, GraphFactory.createPartitionedGraphFromPlipData(chain,
                    plipDocument,
                    inflation,
                    GraphFactory.WeightingScheme.CLASSIFIED));
            inSilicoData.put("ene-" + i, GraphFactory.createPartitionedGraphFromPlipData(chain,
                    plipDocument,
                    inflation,
                    GraphFactory.WeightingScheme.ENERGETIC));
        }

        System.out.println(GraphVisualizer.composeJsonRepresentation(chain, exp, inSilicoData));
    }
}
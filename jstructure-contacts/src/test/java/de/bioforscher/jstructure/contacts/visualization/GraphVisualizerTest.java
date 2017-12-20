package de.bioforscher.jstructure.contacts.visualization;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import de.bioforscher.jstructure.contacts.ProteinGraphFactory;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.impl.MarkovClusteringAlgorithm;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
import org.jsoup.Jsoup;
import org.junit.Ignore;
import org.junit.Test;

import java.io.InputStream;
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

    @Test
    @Ignore
    public void shouldComposeGraphRepresentationFor1brr_A() throws JsonProcessingException {
        Structure structure = StructureParser.source("1brr")
                .minimalParsing(true)
                .parse();
        Chain chain = structure.select()
                .chainName("A")
                .asChain();

        new PLIPIntraMolecularAnnotator().process(structure);
        Graph<AminoAcid> graph = ProteinGraphFactory.createProteinGraph(chain, ProteinGraphFactory.InteractionScheme.SALENTIN2015);

        System.out.println(new ObjectMapper().writeValueAsString(new JsonGraph(graph)));
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

        List<AminoAcid> earlyFoldingResidues = TestUtils.getResourceAsStream("start2fold/" + id + ".residues")
                .filter(line -> !line.startsWith("#"))
                .filter(line -> line.endsWith("EARLY"))
                .map(line -> line.split(";")[0])
                .mapToInt(Integer::valueOf)
                .mapToObj(resNum -> chain.select()
                        .aminoAcids()
                        .residueNumber(resNum)
                        .asAminoAcid())
                .collect(Collectors.toList());

        List<Double> dynamine = TestUtils.getResourceAsStream("dynamine/" + id + ".pred")
                .filter(line -> line.contains("\t"))
                .map(line -> line.split("\t")[1])
                .mapToDouble(Double::valueOf)
                .boxed()
                .collect(Collectors.toList());
        List<Double> efoldmine = TestUtils.getResourceAsStream("efoldmine/" + id + ".pred")
                .mapToDouble(Double::valueOf)
                .boxed()
                .collect(Collectors.toList());

        InputStream expLines = TestUtils.getResourceAsInputStream("exp/" + id + ".exp");

        String plipDocument = TestUtils.getResourceAsStream("plip/" + id + ".plip")
                .collect(Collectors.joining(System.lineSeparator()));
        new PLIPIntraMolecularAnnotator().process(chain, Jsoup.parse(plipDocument));
        Graph<AminoAcid> graph = ProteinGraphFactory.createProteinGraph(chain, ProteinGraphFactory.InteractionScheme.SALENTIN2015);

        PartitionedGraph<AminoAcid> exp = ProteinGraphFactory.Partitioned.fromDefinitionFile(graph, expLines);
        Map<String, PartitionedGraph<AminoAcid>> inSilicoData = new TreeMap<>();

        inSilicoData.put("NetCarto", ProteinGraphFactory.Partitioned.fromNetCartoFile(graph, TestUtils.getResourceAsInputStream("netcarto/" + id + "_plip.modules.dat")));

        for(double inflation = 1.2; inflation <= /*2.01*/ 1.21; inflation = inflation + 0.05) {
            int i = (int) Math.round(inflation * 100);
            System.out.println(i);
            inSilicoData.put("cla-" + i, ProteinGraphFactory.Partitioned.fromMCL(graph, MarkovClusteringAlgorithm.DEFAULT_EXPAND_FACTOR, inflation));
        }

        System.out.println(GraphVisualizer.composeJsonRepresentation(chain,
                earlyFoldingResidues,
                dynamine,
                efoldmine,
                exp,
                inSilicoData));
    }
}
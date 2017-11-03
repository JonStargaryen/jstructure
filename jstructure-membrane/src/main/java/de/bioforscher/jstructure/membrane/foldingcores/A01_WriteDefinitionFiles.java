package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.impl.MarkovClusteringAlgorithm;
import de.bioforscher.jstructure.membrane.GraphFactory;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.visualization.GraphVisualizer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

public class A01_WriteDefinitionFiles {
    public static void main(String[] args) {
        MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("order"))
                .forEach(A01_WriteDefinitionFiles::handleFile);
    }

    private static void handleFile(Path path) {
        System.out.println(path.toFile().getName());
        Path basePath = Paths.get("/home/bittrich/tmp/");
        Path plipPath = path.getParent().getParent().getParent().resolve("modularity").resolve("plip");

        String id = MembraneConstants.lines(path)
                .filter(line -> line.contains("pdb:"))
                .map(line -> line.split(": ")[1])
                .findFirst()
                .get();
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];

        Structure structure = StructureParser.source(pdbId)
                .minimalParsing(true)
                .parse();
        Chain chain = structure.select()
                .chainName(chainId)
                .asChain();

        Path plipFile = plipPath.resolve(id + ".plip");
        String plipDocument = MembraneConstants.lines(plipFile)
                .collect(Collectors.joining(System.lineSeparator()));
        new PLIPIntraMolecularAnnotator().process(chain, Jsoup.parse(plipDocument));
        Graph<AminoAcid> graph = GraphFactory.createProteinGraph(chain, GraphFactory.InteractionScheme.SALENTIN2015);
        PartitionedGraph<AminoAcid> exp = GraphFactory.Partitioned.fromDefinitionFile(graph, path);
        MembraneConstants.write(basePath.resolve(id + "_exp.mod"), GraphVisualizer.composeDefinitionString(exp));

        for(double inflation = 1.2; inflation <= 2.0; inflation = inflation + 0.05) {
            int i = (int) Math.round(inflation * 100);
            System.out.println(i);
            PartitionedGraph<AminoAcid> partitionedGraph = GraphFactory.Partitioned.fromMCL(graph, MarkovClusteringAlgorithm.DEFAULT_EXPAND_FACTOR, inflation);
            MembraneConstants.write(basePath.resolve(id + "_" + i + ".mod"),
                    GraphVisualizer.composeDefinitionString(partitionedGraph));
        }
    }
}

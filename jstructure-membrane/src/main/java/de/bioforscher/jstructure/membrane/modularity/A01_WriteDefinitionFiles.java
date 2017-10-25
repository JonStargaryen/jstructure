package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.modularity.visualization.GraphVisualizer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.nio.file.Path;
import java.nio.file.Paths;

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

        PartitionedGraph<AminoAcid> exp = GraphFactory.createPartitionedGraphFromDefinitionFile(chain, path);
        MembraneConstants.write(basePath.resolve(id + "_exp.mod"), GraphVisualizer.composeDefinitionString(exp));

        for(double inflation = 1.2; inflation <= 2.0; inflation = inflation + 0.05) {
            int i = (int) Math.round(inflation * 100);
            System.out.println(i);
            Path plipFile = plipPath.resolve(id + ".plip");
            PartitionedGraph<AminoAcid> unweightedGraph = GraphFactory.createPartitionedGraphFromPlipData(chain,
                    plipFile,
                    inflation,
                    GraphFactory.WeightingScheme.UNWEIGHTED);
            PartitionedGraph<AminoAcid> classifiedGraph = GraphFactory.createPartitionedGraphFromPlipData(chain,
                    plipFile,
                    inflation,
                    GraphFactory.WeightingScheme.CLASSIFIED);
            PartitionedGraph<AminoAcid> energeticGraph = GraphFactory.createPartitionedGraphFromPlipData(chain,
                    plipFile,
                    inflation,
                    GraphFactory.WeightingScheme.ENERGETIC);
            MembraneConstants.write(basePath.resolve(id + "_unweighted_" + i + ".mod"),
                    GraphVisualizer.composeDefinitionString(unweightedGraph));
            MembraneConstants.write(basePath.resolve(id + "_classified_" + i + ".mod"),
                    GraphVisualizer.composeDefinitionString(classifiedGraph));
            MembraneConstants.write(basePath.resolve(id + "_energetic_" + i + ".mod"),
                    GraphVisualizer.composeDefinitionString(energeticGraph));
        }
    }
}

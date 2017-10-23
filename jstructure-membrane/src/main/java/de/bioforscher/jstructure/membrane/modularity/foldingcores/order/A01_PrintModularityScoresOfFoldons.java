package de.bioforscher.jstructure.membrane.modularity.foldingcores.order;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.PartitioningScorer;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.modularity.GraphFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.nio.file.Path;
import java.util.stream.Collectors;

/**
 *
 */
public class A01_PrintModularityScoresOfFoldons {
    public static void main(String[] args) {
        System.out.println("order;size;intra;inter;score;loop;asa;ep");
        MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("order"))
                .forEach(A01_PrintModularityScoresOfFoldons::handleFile);
    }

    private static void handleFile(Path path) {
        String id = MembraneConstants.lines(path)
                .filter(line -> line.contains("pdb:"))
                .map(line -> line.split(": ")[1])
                .findFirst()
                .get();
        System.out.println(id);
        Structure structure = StructureParser.source(id.split("_")[0])
                .minimalParsing(true)
                .parse();
        Chain chain = structure.select()
                .chainName(id.split("_")[1])
                .asChain();
        new DictionaryOfProteinSecondaryStructure().process(structure);
        new LoopFractionCalculator().process(structure);
        new AccessibleSurfaceAreaCalculator().process(structure);
        new EnergyProfileCalculator().process(structure);
        String plipDocument = MembraneConstants.lines(path.getParent()
                .getParent()
                .getParent()
                .resolve("modularity")
                .resolve("plip")
                .resolve(id + ".plip"))
                .collect(Collectors.joining(System.lineSeparator()));
        Graph<AminoAcid> graph = GraphFactory.createGraphFromPlipDocument(chain, plipDocument, GraphFactory.WeightingScheme.CLASSIFIED);
        PartitionedGraph<AminoAcid> partitionedGraph = GraphFactory.createPartitionedGraphFromDefinitionFile(graph, path);
        partitionedGraph.getModules().stream()
                .map(module -> module.getIdentifier() + ";" +
                        module.getNodes().size() + ";" +
                        partitionedGraph.getEdges().stream()
                                .filter(edge -> module.containsNode(edge.getLeft()) && module.containsNode(edge.getRight()))
                                .count() + ";" +
                        partitionedGraph.getEdges().stream()
                                .filter(edge -> (module.containsNode(edge.getLeft()) && !module.containsNode(edge.getRight())) || (!module.containsNode(edge.getLeft()) && module.containsNode(edge.getRight())))
                                .count()+ ";" +
                        StandardFormat.format(PartitioningScorer.modularityScore(partitionedGraph, module)) + ";" +
                        StandardFormat.format(module.getNodes().stream()
                                .map(aminoAcid -> aminoAcid.getFeature(LoopFraction.class))
                                .mapToDouble(LoopFraction::getLoopFraction)
                                .average()
                                .getAsDouble()) + ";" +
                        StandardFormat.format(module.getNodes().stream()
                                .map(aminoAcid -> aminoAcid.getFeature(AccessibleSurfaceArea.class))
                                .mapToDouble(AccessibleSurfaceArea::getRelativeAccessibleSurfaceArea)
                                .average()
                                .getAsDouble()) + ";" +
                        StandardFormat.format(module.getNodes().stream()
                                .map(aminoAcid -> aminoAcid.getFeature(EnergyProfile.class))
                                .mapToDouble(EnergyProfile::getSolvationEnergy)
                                .average()
                                .getAsDouble()))
                .forEach(System.out::println);
    }
}

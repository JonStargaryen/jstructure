package de.bioforscher.jstructure.mathematics.graph.clustering;

import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.model.SetOperations;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Compares the similarity of two clustering solutions.
 * reference: Wagner and Wagner, 2007 - Comparing Clusterings - An Overview
 */
public class ClusteringScorer {
    /**
     * Related to Meila-Heckerman- and Maximum-Match-Measure.
     * @param clustering1
     * @param clustering2
     * @param <N>
     * @return
     */
    public static <N> double naiveScore(PartitionedGraph<N> clustering1, PartitionedGraph<N> clustering2) {
        // determine more refined clustering
        PartitionedGraph<N> fineClustering = clustering2.getNumberOfModules() > clustering1.getNumberOfModules() ? clustering2 : clustering1;
        PartitionedGraph<N>  coarseClustering = clustering2.getNumberOfModules() > clustering1.getNumberOfModules() ? clustering1 : clustering2;

        // determine total number of elements
        int totalNumberOfElements = extractAllReferencedNodes(fineClustering, coarseClustering).size();

        // for each of fine clusters
        return fineClustering.getModules().stream()
                // determine maximum number of shared elements with a coarse set
                .mapToInt(fineModule -> coarseClustering.getModules().stream()
                        .map(coarseModule -> SetOperations.intersection(fineModule.getNodes(), coarseModule.getNodes()))
                        .mapToInt(Collection::size)
                        .max()
                        .orElse(0))
                .sum() / (double) totalNumberOfElements;
    }

    public static <N> double chiSquaredCoefficient(PartitionedGraph<N> clustering1, PartitionedGraph<N> clustering2) {
        double chi = 0;
        int totalNumberOfElements = extractAllReferencedNodes(clustering1, clustering2).size();

        for(int i = 0; i < clustering1.getNumberOfModules(); i++) {
            Module<N> module1 = clustering1.getModules().get(i);
            int size1 = module1.getNumberOfNodes();

            for(int j = 0; j < clustering2.getNumberOfModules(); j++) {
                Module<N> module2 = clustering2.getModules().get(j);
                int size2 = module2.getNumberOfNodes();

                double eij = (size1 * size2) / (double) totalNumberOfElements;
                double mij = SetOperations.union(module1.getNodes(), module2.getNodes()).size();
                double me = mij - eij;

                chi += (me * me) / eij;
            }
        }

        return chi;
    }

    public static <N> double randIndex(PartitionedGraph<N> clustering1, PartitionedGraph<N> clustering2) {
        List<N> nodes = extractAllReferencedNodes(clustering1, clustering2);
        int totalNumberOfElements = nodes.size();

        int n11 = 0; // in same cluster in both solutions
        int n00 = 0; // in different clusters in both solutions

        for(int i = 0; i < nodes.size() - 1; i++) {
            N node1 = nodes.get(i);
            Module module11 = clustering1.getModuleOf(node1);
            Module module12 = clustering2.getModuleOf(node1);
            for(int j = i + 1; j < nodes.size(); j++) {
                N node2 = nodes.get(j);
                Module module21 = clustering1.getModuleOf(node2);
                Module module22 = clustering2.getModuleOf(node2);

                // same module in both solutions
                if(module11.equals(module21) && module21.equals(module22)) {
                    n11++;
                } else if(!module11.equals(module21) && !module12.equals(module22)) {
                    n00++;
                }
            }
        }

        return 2 * (n11 + n00) / (double) (totalNumberOfElements * (totalNumberOfElements - 1));
    }

    private static <N> List<N> extractAllReferencedNodes(PartitionedGraph<N> clustering1, PartitionedGraph<N> clustering2) {
        return Stream.of(clustering1, clustering2)
                .map(PartitionedGraph::getNodes)
                .flatMap(Collection::stream)
                .distinct()
                .collect(Collectors.toList());
    }
}

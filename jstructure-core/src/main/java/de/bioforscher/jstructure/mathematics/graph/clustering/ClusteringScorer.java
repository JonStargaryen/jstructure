package de.bioforscher.jstructure.mathematics.graph.clustering;

import de.bioforscher.jstructure.model.SetOperations;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Compares the similarity of two clustering solutions.
 * reference: Wagner and Wagner, 2007 - Comparing Clusterings - An Overview
 */
public class ClusteringScorer {
    public static <N> double naiveScore(List<Module<N>> clustering1, List<Module<N>> clustering2) {
        // determine more refined clustering
        List<Module<N>> fineClustering = clustering2.size() > clustering1.size() ? clustering2 : clustering1;
        List<Module<N>> coarseClustering = clustering2.size() > clustering1.size() ? clustering1 : clustering2;

        // determine total number of elements
        int totalNumberOfElements = extractNodes(fineClustering, coarseClustering).size();

        // for each of fine clusters
        return fineClustering.stream()
                // determine maximum number of shared elements with a coarse set
                .mapToInt(fineModule -> coarseClustering.stream()
                        .map(coarseModule -> SetOperations.intersection(fineModule.getNodes(), coarseModule.getNodes()))
                        .mapToInt(Collection::size)
                        .max()
                        .orElse(0))
                .sum() / (double) totalNumberOfElements;
    }

    public static <N> double chiSquaredCoefficient(List<Module<N>> clustering1, List<Module<N>> clustering2) {
        double chi = 0;
        int totalNumberOfElements = extractNodes(clustering1, clustering2).size();

        for(int i = 0; i < clustering1.size(); i++) {
            Module<N> module1 = clustering1.get(i);
            int size1 = module1.getNumberOfNodes();

            for(int j = 0; j < clustering2.size(); j++) {
                Module<N> module2 = clustering2.get(j);
                int size2 = module2.getNumberOfNodes();

                double eij = (size1 * size2) / (double) totalNumberOfElements;
                double mij = SetOperations.union(module1.getNodes(), module2.getNodes()).size();
                double me = mij - eij;

                chi += (me * me) / eij;
            }
        }

        return chi;
    }

    public static <N> double randIndex(List<Module<N>> clustering1, List<Module<N>> clustering2) {
        List<N> nodes = extractNodes(clustering1, clustering2);
        int totalNumberOfElements = nodes.size();

        int n11 = 0; // in same cluster in both solutions
        int n00 = 0; // in different clusters in both solutions

        for(int i = 0; i < nodes.size() - 1; i++) {
            N node1 = nodes.get(i);
            Optional<Module<N>> module11 = findModuleOfNode(node1, clustering1);
            Optional<Module<N>> module12 = findModuleOfNode(node1, clustering2);
            for(int j = i + 1; j < nodes.size(); j++) {
                N node2 = nodes.get(j);
                Optional<Module<N>> module21 = findModuleOfNode(node2, clustering1);
                Optional<Module<N>> module22 = findModuleOfNode(node2, clustering2);

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

    private static <N> Optional<Module<N>> findModuleOfNode(N node, List<Module<N>> clustering) {
        return clustering.stream()
                .filter(module -> module.containsNode(node))
                .findFirst();
    }

    private static <N> List<N> extractNodes(List<Module<N>> clustering1, List<Module<N>> clustering2) {
        return Stream.of(clustering1, clustering2)
                .flatMap(Collection::stream)
                .map(Module::getNodes)
                .flatMap(Collection::stream)
                .distinct()
                .collect(Collectors.toList());
    }
}

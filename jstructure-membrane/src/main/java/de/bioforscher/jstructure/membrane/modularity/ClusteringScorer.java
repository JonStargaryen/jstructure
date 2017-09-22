package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.model.SetOperations;

import java.util.Collection;
import java.util.List;

/**
 * Compares the similarity of two clustering solutions.
 */
public class ClusteringScorer {
    public static double naiveScore(List<Module> clustering1, List<Module> clustering2) {
        // determine more refined clustering
        List<Module> fineClustering = clustering2.size() > clustering1.size() ? clustering2 : clustering1;
        List<Module> coarseClustering = clustering2.size() > clustering1.size() ? clustering1 : clustering2;

        // determine total number of elements
        int totalNumberOfElements = fineClustering.stream()
                .mapToInt(Module::getSize)
                .sum();

        // for each of fine clusters
        return fineClustering.stream()
                // determine maximum number of shared elements with a coarse set
                .mapToInt(fineModule -> coarseClustering.stream()
                        .map(coarseModule -> SetOperations.intersection(fineModule.getNodeNames(), coarseModule.getNodeNames()))
                        .mapToInt(Collection::size)
                        .max()
                        .orElse(0))
                .sum() / (double) totalNumberOfElements;
    }

    /**
     * reference: Wagner, 2007 - Comparing Clusterings - An Overview
     * @param clustering1
     * @param clustering2
     * @return
     */
    public static double chiSquaredCoefficient(List<Module> clustering1, List<Module> clustering2) {
        double chi = 0;
        int totalNumberOfElements = clustering1.stream()
                .mapToInt(Module::getSize)
                .sum();

        for(int i = 0; i < clustering1.size(); i++) {
            Module module1 = clustering1.get(i);
            int size1 = module1.getSize();

            for(int j = 0; j < clustering2.size(); j++) {
                Module module2 = clustering2.get(j);
                int size2 = module2.getSize();

                double eij = (size1 * size2) / (double) totalNumberOfElements;
                double mij = SetOperations.union(module1.getNodeNames(), module2.getNodeNames()).size();
                double me = mij - eij;

                chi += (me * me) / eij;
            }
        }

        return chi;
    }
}

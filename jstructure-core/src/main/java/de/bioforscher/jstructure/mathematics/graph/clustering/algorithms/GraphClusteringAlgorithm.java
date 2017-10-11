package de.bioforscher.jstructure.mathematics.graph.clustering.algorithms;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.clustering.Module;

import java.util.List;

public interface GraphClusteringAlgorithm {
    <E> List<Module<E>> clusterGraph(Graph<E> graph);
}

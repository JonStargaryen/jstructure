package de.bioforscher.jstructure.mathematics.graph.clustering.algorithms;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;

public interface GraphClusteringAlgorithm {
    <N> PartitionedGraph<N> clusterGraph(Graph<N> graph);
}

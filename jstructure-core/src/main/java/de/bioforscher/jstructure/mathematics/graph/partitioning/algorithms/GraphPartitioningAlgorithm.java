package de.bioforscher.jstructure.mathematics.graph.partitioning.algorithms;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;

public interface GraphPartitioningAlgorithm {
    <N> PartitionedGraph<N> partitionGraph(Graph<N> graph);
}

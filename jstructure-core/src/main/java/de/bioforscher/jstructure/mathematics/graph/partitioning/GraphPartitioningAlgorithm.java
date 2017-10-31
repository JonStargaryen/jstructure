package de.bioforscher.jstructure.mathematics.graph.partitioning;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;

/**
 * Algorithms used to partition graphs into a set of {@link Module} instances.
 */
public interface GraphPartitioningAlgorithm {
    /**
     * Partitions a graph into disjoint sets of subgraphs.
     * @param graph the input graph
     * @param <N> the node type
     * @return a partitioned graph (optimal in the eyes of the employed algorithm)
     */
    <N> PartitionedGraph<N> partitionGraph(Graph<N> graph);
}

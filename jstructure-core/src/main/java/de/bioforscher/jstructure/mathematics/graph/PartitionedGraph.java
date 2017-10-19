package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;
import de.bioforscher.jstructure.mathematics.graph.partitioning.PartitioningScorer;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class PartitionedGraph<N> extends Graph<N> {
    private final List<Module<N>> modules;

    public PartitionedGraph(Graph<N> graph, List<Module<N>> modules) {
        super(graph);
        this.modules = modules;
    }

    public Optional<Module<N>> getModuleOf(N node) {
        return modules.stream()
                .filter(module -> module.getNodes().contains(node))
                .findFirst();
    }

    public List<Module<N>> getModules() {
        return modules;
    }

    public int getNumberOfModules() {
        return modules.size();
    }

    public List<String> getModuleIdentifiers() {
        return modules.stream()
                .map(Module::getIdentifier)
                .collect(Collectors.toList());
    }

    public double getModularityScore() {
        return PartitioningScorer.modularityScore(this);
    }
}

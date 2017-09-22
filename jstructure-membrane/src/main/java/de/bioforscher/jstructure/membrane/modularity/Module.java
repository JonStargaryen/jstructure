package de.bioforscher.jstructure.membrane.modularity;

import java.util.List;

/**
 * Represents a module or cluster or partition or community in a graph or network. Basically just a association of
 * amino acids and the contract, that for a chain each amino acid can only be assigned to one module.
 */
public class Module {
    private final List<String> nodeNames;

    public Module(List<String> nodeNames) {
        this.nodeNames = nodeNames;
    }

    public List<String> getNodeNames() {
        return nodeNames;
    }

    public int getSize() {
        return nodeNames.size();
    }
}

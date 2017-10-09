package de.bioforscher.jstructure.membrane.modularity;

import java.util.List;

/**
 * Represents a module or cluster or partition or community in a graph or network. Basically just a association of
 * amino acids and the contract, that for a chain each amino acid can only be assigned to one module.
 */
public class Module {
    private final String id;
    private final List<String> nodeNames;

    public Module(String id, List<String> nodeNames) {
        this.id = id;
        this.nodeNames = nodeNames;
    }

    public String getId() {
        return id;
    }

    public List<String> getNodeNames() {
        return nodeNames;
    }

    public int getSize() {
        return nodeNames.size();
    }
}

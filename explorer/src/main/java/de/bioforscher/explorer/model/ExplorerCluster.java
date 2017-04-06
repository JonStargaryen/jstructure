package de.bioforscher.explorer.model;

import java.util.List;

/**
 * Created by bittrich on 4/6/17.
 */
public class ExplorerCluster {
    private List<ExplorerChain> explorerChains;
    private ExplorerAlignment explorerAlignment;

    public ExplorerCluster() {

    }

    public ExplorerCluster(List<ExplorerChain> explorerChains, ExplorerAlignment alignment) {
        this.explorerChains = explorerChains;
        this.explorerAlignment = alignment;
    }

    public List<ExplorerChain> getExplorerChains() {
        return explorerChains;
    }

    public ExplorerAlignment getExplorerAlignment() {
        return explorerAlignment;
    }
}

package de.bioforscher.explorer.membrane.model;

import java.util.List;

/**
 * Container of stuff to persist.
 * Created by bittrich on 3/23/17.
 */
public class ExplorerCluster {
    private ExplorerAlignment alignment;
    private List<ExplorerChain> chains;

    public ExplorerCluster(ExplorerAlignment alignment, List<ExplorerChain> chains) {
        this.alignment = alignment;
        this.chains = chains;
    }

    public ExplorerAlignment getAlignment() {
        return alignment;
    }

    public List<ExplorerChain> getChains() {
        return chains;
    }
}

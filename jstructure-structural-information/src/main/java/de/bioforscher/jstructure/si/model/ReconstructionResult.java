package de.bioforscher.jstructure.si.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;

public class ReconstructionResult {
    private final List<Chain> chains;
    private final long time;

    public ReconstructionResult(List<Chain> chains, long time) {
        this.chains = chains;
        this.time = time;
    }

    public List<Chain> getChains() {
        return chains;
    }

    public long getTime() {
        return time;
    }
}

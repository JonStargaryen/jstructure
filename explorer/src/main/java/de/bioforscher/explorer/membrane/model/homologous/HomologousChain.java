package de.bioforscher.explorer.membrane.model.homologous;

/**
 * Represents an aligned chain's sequence.
 * Created by bittrich on 3/16/17.
 */
public class HomologousChain {
    private String id, sequence, comp;

    public HomologousChain() {
    }

    public HomologousChain(String id, String sequence) {
        this.id = id;
        this.sequence = sequence;
    }

    public String getId() {
        return id;
    }

    public String getSequence() {
        return sequence;
    }

    @Override
    public String toString() {
        return id + System.lineSeparator() + sequence;
    }
}

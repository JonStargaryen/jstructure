package de.bioforscher.explorer.membrane.model;

/**
 * An aligned sequence.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerSequence {
    private String id, sequence;

    public ExplorerSequence() {
    }

    public ExplorerSequence(String id, String sequence) {
        this.id = id;
        this.sequence = sequence;
    }

    public String getId() {
        return id;
    }

    public String getSequence() {
        return sequence;
    }

    public void setId(String id) {
        this.id = id;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
}

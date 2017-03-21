package de.bioforscher.explorer.membrane.model;

/**
 * No real {@link ExplorerGroup}, no real <code>char</code>.
 * Created by bittrich on 3/21/17.
 */
public class SequencePosition {
    private int resn;
    private String olc;

    public SequencePosition() {
    }

    public SequencePosition(int resn, String olc) {
        this.resn = resn;
        this.olc = olc;
    }

    public int getResn() {
        return resn;
    }

    public String getOlc() {
        return olc;
    }
}

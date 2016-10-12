package de.bioforscher.jstructure.feature.sse;

/**
 * A Ladder is a set of one or more consecutive bridges of identical type. A
 * Bridge is a Ladder of length one.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class Ladder {
    private BridgeType btype;
    // another ladder with lower index connected to this
    private int connectedFrom;
    // another ladder with higher index connected to this
    private int connectedTo;
    private int from; // start of the first strand

    private int lfrom; // start of the second strand

    private int lto; // end of the second strand

    private int to; // end of the first strand

    public BridgeType getBtype() {
        return this.btype;
    }

    public int getConnectedFrom() {
        return this.connectedFrom;
    }

    public int getConnectedTo() {
        return this.connectedTo;
    }

    public int getFrom() {
        return this.from;
    }

    public int getLfrom() {
        return this.lfrom;
    }

    public int getLto() {
        return this.lto;
    }

    public int getTo() {
        return this.to;
    }

    public void setBtype(BridgeType btype) {
        this.btype = btype;
    }

    public void setConnectedFrom(int connectedFrom) {
        this.connectedFrom = connectedFrom;
    }

    public void setConnectedTo(int connectedTo) {
        this.connectedTo = connectedTo;
    }

    public void setFrom(int from) {
        this.from = from;
    }

    public void setLfrom(int lfrom) {
        this.lfrom = lfrom;
    }

    public void setLto(int lto) {
        this.lto = lto;
    }

    public void setTo(int to) {
        this.to = to;
    }

    @Override
    public String toString() {
        return "Ladder [from=" + from + ", to=" + to + ", lfrom=" + lfrom + ", lto=" + lto + ", btype=" + btype + ", connectedTo=" + connectedTo + ", connectedFrom=" + connectedFrom + "]";
    }
}
package de.bioforscher.jstructure.membrane.modularity.visualization;

public class JsonRange {
    private final int x;
    private final int y;

    public JsonRange(int start, int end) {
        this.x = start;
        this.y = end;
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }
}

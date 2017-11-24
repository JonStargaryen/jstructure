package de.bioforscher.jstructure.contacts.visualization;

import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class JsonEdge {
    private final int source;
    private final int target;
    private final double value;

    public JsonEdge(Edge<AminoAcid> edge) {
        this.source = edge.getLeft().getResidueIndex();
        this.target = edge.getRight().getResidueIndex();
        this.value = edge.getWeight();
    }

    public int getSource() {
        return source;
    }

    public int getTarget() {
        return target;
    }

    public double getValue() {
        return value;
    }
}

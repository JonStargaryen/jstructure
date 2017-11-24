package de.bioforscher.jstructure.contacts.visualization;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.stream.Collectors;

public class JsonGraph {
    private final List<JsonNode> nodes;
    private final List<JsonEdge> edges;

    public JsonGraph(Graph<AminoAcid> graph) {
        this.nodes = graph.getNodes()
                .stream()
                .map(JsonNode::new)
                .collect(Collectors.toList());
        this.edges = graph.getEdges()
                .stream()
                .map(JsonEdge::new)
                .collect(Collectors.toList());
    }

    public List<JsonNode> getNodes() {
        return nodes;
    }

    public List<JsonEdge> getEdges() {
        return edges;
    }
}

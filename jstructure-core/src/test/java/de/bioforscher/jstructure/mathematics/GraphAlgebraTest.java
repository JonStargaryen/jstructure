package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * http://www.vogella.com/tutorials/JavaAlgorithmsDijkstra/article.html
 */
public class GraphAlgebraTest {
    private List<String> nodes;
    private List<Edge<String>> edges;
    private Graph<String> graph;

    @Before
    public void setup() {
        nodes = new ArrayList<>();
        edges = new ArrayList<>();
        for (int i = 0; i < 11; i++) {
            String location = "Node_" + i;
            nodes.add(location);
        }

        addLane("Edge_0", 0, 1, 85);
        addLane("Edge_1", 0, 2, 217);
        addLane("Edge_2", 0, 4, 173);
        addLane("Edge_3", 2, 6, 186);
        addLane("Edge_4", 2, 7, 103);
        addLane("Edge_5", 3, 7, 183);
        addLane("Edge_6", 5, 8, 250);
        addLane("Edge_7", 8, 9, 84);
        addLane("Edge_8", 7, 9, 167);
        addLane("Edge_9", 4, 9, 502);
        addLane("Edge_10", 9, 10, 40);
        addLane("Edge_11", 1, 10, 600);

        // Lets check from location Loc_1 to Loc_10
        graph = new Graph<>(nodes, edges);
    }

    @Test
    public void shouldDetermineAveragePathLength() {
        System.out.println(graph.calculate().averageGraphPathLength());
    }

    private void addLane(String laneId,
                         int sourceLocNo,
                         int destLocNo,
                         int duration) {
        Edge<String> lane = new Edge<>(nodes.get(sourceLocNo),
                nodes.get(destLocNo),
                duration);
        edges.add(lane);
    }
}

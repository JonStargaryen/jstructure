package de.bioforscher.jstructure.mathematics.graph.partitioning.algorithms;

import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.stream.Collectors;
import java.util.stream.Stream;

public class MCLTest {
    private MCL mcl;
    private Graph<String> graph;
    private String cat;
    private String hat;
    private String bat;
    private String bit;
    private String fit;
    private String hit;

    @Before
    public void setup() {
        mcl = new MCL();

        cat = "cat";
        hat = "hat";
        bat = "bat";
        bit = "bit";
        fit = "fit";
        hit = "hit";

        Edge<String> cathat = new Edge<>(cat, hat, 0.2);
        Edge<String> hatbat = new Edge<>(hat, bat, 0.16);
        Edge<String> batcat = new Edge<>(bat, cat, 1.0);
        Edge<String> batbit = new Edge<>(bat, bit, 0.125);
        Edge<String> bitfit = new Edge<>(bit, fit, 0.25);
        Edge<String> fithit = new Edge<>(fit, hit, 0.5);
        Edge<String> hitbit = new Edge<>(hit, bit, 0.16);

        graph = new Graph<>(Stream.of(cat,
                hat,
                bat,
                bit,
                fit,
                hit).collect(Collectors.toList()),
            Stream.of(cathat,
                hatbat,
                batcat,
                batbit,
                bitfit,
                fithit,
                hitbit).collect(Collectors.toList()));
    }

    @Test
    public void shouldReportIdentity() {
        Assert.assertTrue(mcl.numberEqual(1, 1, 4));
        Assert.assertTrue(mcl.numberEqual(0.999999, 1, 4));
        Assert.assertTrue(mcl.numberEqual(0.99999, 1, 4));
    }

    @Test
    public void shouldReportNonIdentity() {
        Assert.assertFalse(mcl.numberEqual(1, 2, 4));
        Assert.assertFalse(mcl.numberEqual(1, 1.001, 4));
        Assert.assertFalse(mcl.numberEqual(0.99, 1, 4));
        Assert.assertFalse(mcl.numberEqual(0.999, 1, 4));
    }

    @Test
    public void shouldRunMCL() {
        MCL mcl = new MCL(MCL.DEFAULT_EXPAND_FACTOR,
                1.5,
                MCL.DEFAULT_MULT_FACTOR,
                MCL.DEFAULT_MAX_ITERATIONS,
                Edge::getWeight);

        PartitionedGraph<String> clustering = mcl.partitionGraph(graph);

        Assert.assertEquals("cluster size does not match expectation",
                2,
                clustering.getNumberOfModules());
        Assert.assertTrue("partitioning compromised", clustering.getModules().get(0).containsNode(cat));
        Assert.assertTrue("partitioning compromised", clustering.getModules().get(0).containsNode(hat));
        Assert.assertTrue("partitioning compromised", clustering.getModules().get(0).containsNode(bat));
        Assert.assertTrue("partitioning compromised", clustering.getModules().get(1).containsNode(bit));
        Assert.assertTrue("partitioning compromised", clustering.getModules().get(1).containsNode(hit));
        Assert.assertTrue("partitioning compromised", clustering.getModules().get(1).containsNode(fit));
    }
}
package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.feature.topology.Membrane;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Ensures correct calculations of the distance of a point to a membrane.
 * Created by S on 10.11.2016.
 */
public class MembraneUnitTest {
    private Membrane membrane;
    private double[] pointInMembrane;
    private double[] pointAboveMembrane;
    private double[] pointBelowMembrane;

    @Before
    public void setup() {
        pointAboveMembrane = new double[] { 1, 0, 0 };
        pointBelowMembrane = new double[] { -1, 0, 0 };
        membrane = new Membrane(pointAboveMembrane, pointBelowMembrane, 0, 0);
        pointInMembrane = new double[3];

    }

    @Test
    public void testDistanceToMembraneCenter() {
        // point in plane
        double distanceInPlane = membrane.distanceToMembraneCenter(pointInMembrane);
        Assert.assertEquals(0.0, distanceInPlane, 0.001);

        // point 1 length unit above plane
        double distanceAbovePlane = membrane.distanceToMembraneCenter(pointAboveMembrane);
        Assert.assertEquals(1.0, distanceAbovePlane, 0.001);

        // point 1 length unit below plane
        double distanceBelowPlane = membrane.distanceToMembraneCenter(pointBelowMembrane);
        Assert.assertEquals(-1.0, distanceBelowPlane, 0.001);
    }
}

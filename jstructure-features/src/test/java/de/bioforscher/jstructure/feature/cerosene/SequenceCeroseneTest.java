package de.bioforscher.jstructure.feature.cerosene;

import org.junit.Assert;
import org.junit.Test;

/**
 * Functional and unit tests for {@link SequenceCerosene}.
 * Created by bittrich on 3/14/17.
 */
public class SequenceCeroseneTest {
    @Test
    public void convertBlackToHsv() {
        double[] hsv = SequenceCerosene.rgb2hsv(new double[] { 0, 0, 0 });
        Assert.assertEquals(0.0, hsv[2], 0.0);
    }

    @Test
    public void convertWhiteToHsv() {
        double[] hsv = SequenceCerosene.rgb2hsv(new double[] { 255, 255, 255 });
        Assert.assertEquals(0.0, hsv[1], 0.0);
    }

    @Test
    public void convertGreenToHsv() {
        double[] hsv = SequenceCerosene.rgb2hsv(new double[] { 0, 255, 0 });
        Assert.assertEquals(0.3333, hsv[0], 0.001);
    }

    @Test
    public void convertRedToHsv() {
        double[] hsv = SequenceCerosene.rgb2hsv(new double[] { 255, 0, 0 });
        Assert.assertEquals(0.0, hsv[0], 0.001);
    }
}
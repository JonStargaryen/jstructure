package de.bioforscher.jstructure.mathematics;

import org.junit.Assert;
import org.junit.Test;

/**
 * Test for shared functions by the linear calculate implementation.
 * Created by bittrich on 5/23/17.
 */
public class LinearAlgebraTest {
    @Test
    public void capToInterval() throws Exception {
        // [low: -1, value: 0, high: 1] => 0
        Assert.assertEquals(0.0, LinearAlgebra.capToInterval(-1, 0, 1), 0.0);
        // [low: -1, value: -1, high: 1] => -1
        Assert.assertEquals(-1.0, LinearAlgebra.capToInterval(-1, -1, 1), 0.0);
        // [low: -1, value: -2, high: 1] => -1
        Assert.assertEquals(-1.0, LinearAlgebra.capToInterval(-1, -2, 1), 0.0);
        // [low: -1, value: 2, high: 1] => 1
        Assert.assertEquals(1.0, LinearAlgebra.capToInterval(-1, 2, 1), 0.0);
    }
}
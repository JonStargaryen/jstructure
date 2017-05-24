package de.bioforscher.jstructure.mathematics;

import org.junit.Assert;
import org.junit.Test;

import static util.TestUtils.TOLERANT_ERROR_MARGIN;
import static util.TestUtils.ZERO_VECTOR;

/**
 * Unit tests for linear calculate on double arrays.
 * Created by bittrich on 5/23/17.
 */
public class PrimitiveDoubleArrayLinearAlgebraTest {
    private static final double[] vector1 = new double[] { 1, 2, 3 };
    private static final double[] vector2 = new double[] { 2, 3, 1 };
    private static final double[] orthogonalVector1 = new double[] { 0, 1, 0 };
    private static final double[] orthogonalVector2 = new double[] { 1, 0, 0 };
    private static final double scalar1 = 0.5;
    private static final double scalar2 = 2.0;

    @Test
    public void addAndMultiply() throws Exception {
        Assert.assertArrayEquals(LinearAlgebra.on(vector1).add(vector1).getValue(),
            LinearAlgebra.on(vector1).multiply(scalar2).getValue(),
            TOLERANT_ERROR_MARGIN);
        Assert.assertArrayEquals(LinearAlgebra.on(vector2).add(vector2).getValue(),
                LinearAlgebra.on(vector2).multiply(scalar2).getValue(),
                TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void angle() throws Exception {
        // vector is the same => 0°
        Assert.assertEquals(0.0, LinearAlgebra.on(vector1).angle(vector1), 0.0);
        // vector is orthogonal => 90°
        Assert.assertEquals(90.0, LinearAlgebra.on(orthogonalVector1).angle(orthogonalVector2), 0.0);
        // realistic case
        Assert.assertEquals(Math.toDegrees(0.666946), LinearAlgebra.on(vector1).angle(vector2), TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void distance() throws Exception {
        // zero distance for same point
        Assert.assertEquals(0.0, LinearAlgebra.on(vector1).distance(vector1), 0.0);
        // should report norm
        Assert.assertEquals(LinearAlgebra.on(vector1).norm(),
                LinearAlgebra.on(vector1).distance(ZERO_VECTOR), 0.0);
    }

    @Test
    public void divide() throws Exception {
        Assert.assertArrayEquals(LinearAlgebra.on(vector1).multiply(scalar1).getValue(),
                LinearAlgebra.on(vector1).divide(scalar2).getValue(),
                0.0);
    }

    @Test
    public void dotProduct() throws Exception {
        // orthogonal vectors => 0.0
        Assert.assertEquals(0.0, LinearAlgebra.on(orthogonalVector1).dotProduct(orthogonalVector2), 0.0);
    }

    @Test
    public void normalize() throws Exception {
        Assert.assertEquals(1.0, LinearAlgebra.on(vector1).normalize().norm(), 0.0);
    }

    @Test
    public void vectorProduct() throws Exception {
        //TODO impl
    }
}
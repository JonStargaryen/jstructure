package unit.mathematics;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by S on 24.10.2016.
 */
public class LinearAlgebra3DUnitTest {
    private static final double[] vector1 = new double[] { 1, 2, 3 };
    private static final double[] vector2 = new double[] { 2, 3, 1 };
    private static final double[] orthogonalVector1 = new double[] { 0, 1, 0 };
    private static final double[] orthogonalVector2 = new double[] { 1, 0, 0 };
    private static final double scalar1 = 0.5;
    private static final double scalar2 = 2.0;

    private static final double[] ZERO_VECTOR = new double[] { 0, 0, 0 };
    public static final double TOLERANT_ERROR_MARGIN = 0.001;

    @Test
    public void testAdditionAndMultiplication() {
        Assert.assertArrayEquals(LinearAlgebra3D.add(vector1, vector1),
                LinearAlgebra3D.multiply(vector1, scalar2),
                TOLERANT_ERROR_MARGIN);
        Assert.assertArrayEquals(LinearAlgebra3D.add(vector2, vector2),
                LinearAlgebra3D.multiply(vector2, scalar2),
                TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void testAngle() {
        // vector is the same => 0°
        Assert.assertEquals(0.0, LinearAlgebra3D.angle(vector1, vector1), 0.0);
        // vector is orthogonal => 90°
        Assert.assertEquals(90.0, LinearAlgebra3D.angle(orthogonalVector1, orthogonalVector2), 0.0);
        // realistic case
        Assert.assertEquals(Math.toDegrees(0.666946), LinearAlgebra3D.angle(vector1, vector2), TOLERANT_ERROR_MARGIN);
    }

    @Test
    public void testCapToInterval() {
        // [low: -1, value: 0, high: 1] => 0
        Assert.assertEquals(0.0, LinearAlgebra3D.capToInterval(-1, 0, 1), 0.0);
        // [low: -1, value: -1, high: 1] => -1
        Assert.assertEquals(-1.0, LinearAlgebra3D.capToInterval(-1, -1, 1), 0.0);
        // [low: -1, value: -2, high: 1] => -1
        Assert.assertEquals(-1.0, LinearAlgebra3D.capToInterval(-1, -2, 1), 0.0);
        // [low: -1, value: 2, high: 1] => 1
        Assert.assertEquals(1.0, LinearAlgebra3D.capToInterval(-1, 2, 1), 0.0);
    }

    @Test
    public void testDistance() {
        // zero distance for same point
        Assert.assertEquals(0.0, LinearAlgebra3D.distance(vector1, vector1), 0.0);
        // should report norm
        Assert.assertEquals(LinearAlgebra3D.norm(vector1),
                LinearAlgebra3D.distance(vector1, ZERO_VECTOR), 0.0);
    }

    @Test
    public void testDivide() {
        Assert.assertArrayEquals(LinearAlgebra3D.multiply(vector1, 0.5), LinearAlgebra3D.divide(vector1, 2), 0.0);
    }

    @Test
    public void testDotProduct() {
        // orthogonal vectors => 0.0
        Assert.assertEquals(0.0, LinearAlgebra3D.dotProduct(orthogonalVector1, orthogonalVector2), 0.0);
    }

    @Test
    public void testNormalize() {
        Assert.assertEquals(1.0, LinearAlgebra3D.norm(LinearAlgebra3D.normalize(vector1)), 0.0);
    }

    @Test
    public void testTorsionAngle() {
        //TODO implement
    }

    @Test
    public void testVectorProduct() {
        //TODO implement
    }
}
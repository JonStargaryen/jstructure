package de.bioforscher.jstructure.mathematics;

/**
 * Provides access to common algebraic computations for 3D vectors. No function guarantees that it function properly
 * for other dimensions, nor will any function check for matching requirements etc. By convention none of the provided
 * functions will happen in place (and manipulate the input data), but rather will create a new primitive instance to
 * return the result.
 * @see CoordinateUtils
 * Created by S on 28.09.2016.
 */
public class LinearAlgebra3D {
    private LinearAlgebra3D() {
        // deny instantiation
    }

    /**
     * Add 2 3D vectors and returns their sum.
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the vector sum
     */
    public static double[] add(double[] v1, double[] v2) {
        return new double[] {
            v1[0] + v2[0],
            v1[1] + v2[1],
            v1[2] + v2[2]
        };
    }

    /**
     * Returns the angle in radians between this vector and the vector
     * parameter; the return value is constrained to the range [0,{@link Math#PI}].
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the angle in degrees in the range [-180,180]
     */
    public static double angle(double[] v1, double[] v2) {
        double vDot = dotProduct(v1, v2) / (norm(v1) * norm(v2));
        return Math.toDegrees(Math.acos(capToInterval(-1, vDot, 1)));
    }

    /**
     * Caps a number to a defined interval.
     * @param lowerBound the lowest accepted value x_l
     * @param value the raw value
     * @param upperBound the highest accepted value x_u
     * @return the original value if it lies in the interval [x_l,x_u], else x_l or x_u
     */
    public static double capToInterval(double lowerBound, double value, double upperBound) {
        return value < lowerBound ? lowerBound : value > upperBound ? upperBound : value;
    }

    /**
     * Caps a number to a defined interval.
     * @param lowerBound the lowest accepted value x_l
     * @param value the raw value
     * @param upperBound the highest accepted value x_u
     * @return the original value if it lies in the interval [x_l,x_u], else x_l or x_u
     */
    public static int capToInterval(int lowerBound, int value, int upperBound) {
        return value < lowerBound ? lowerBound : value > upperBound ? upperBound : value;
    }

    /**
     * Calculates the distance between 2 3D points. Calls {@link LinearAlgebra3D#distanceFast(double[], double[])} and
     * takes the square root of the result. You are encouraged to use
     * {@link LinearAlgebra3D#distanceFast(double[], double[])} where appropriate (e.g. to sort distances or to check
     * some distance relative to some threshold: in that case square the threshold value once).
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the euclidean distance
     */
    public static double distance(double[] v1, double[] v2) {
        return Math.sqrt(distanceFast(v1, v2));
    }

    /**
     * Calculates the squared distance between 2 3D points.
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the squared distance
     * @see LinearAlgebra3D#distance(double[], double[])
     */
    public static double distanceFast(double[] v1, double[] v2) {
        return (v1[0] - v2[0])*(v1[0] - v2[0]) + (v1[1] - v2[1])*(v1[1] - v2[1]) + (v1[2] - v2[2])*(v1[2] - v2[2]);
    }

    /**
     * Divides a 3D vector by a scalar.
     * @param v a 3D vector
     * @param scalar a scalar to divide the vector
     * @return the divided vector
     */
    public static double[] divide(double[] v, double scalar) {
        return new double[] {
            v[0] / scalar,
            v[1] / scalar,
            v[2] / scalar
        };
    }

    /**
     * Calculates the scalar product (inner product) of 2 3D vectors.
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the scalar product of both vectors
     */
    public static double dotProduct(double[] v1, double[] v2) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    /**
     * Multiplies each element of the given 3D vector with a scalar.
     * @param v a 3D vector
     * @param scalar a scalar to multiply the vector with
     * @return the multiplied vector
     */
    public static double[] multiply(double[] v, double scalar) {
        return new double[] {
            v[0] * scalar,
            v[1] * scalar,
            v[2] * scalar
        };
    }

    public static double[] multiply(double[] v, double[][] matrix) {
        return new double[] {
            (matrix[0][0] * v[0] + matrix[1][0] * v[1] + matrix[2][0] * v[2]),
            (matrix[0][1] * v[0] + matrix[1][1] * v[1] + matrix[2][1] * v[2]),
            (matrix[0][2] * v[0] + matrix[1][2] * v[1] + matrix[2][2] * v[2])
        };
    }

    /**
     * Calculates the length (the norm) of this vector.
     * @param v a 3D vector
     * @return the norm of that vector
     */
    public static double norm(double[] v) {
        return Math.sqrt(dotProduct(v, v));
    }

    /**
     * Normalizes 3D vectors (i.e. dividing each component by the {@link LinearAlgebra3D#norm(double[])} of the vector).
     * @param v a 3D vector
     * @return the normalized vector
     */
    public static double[] normalize(double[] v) {
        return divide(v, norm(v));
    }

    /**
     * Subtracts 2 3D vectors.
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the difference vector
     */
    public static double[] subtract(double[] v1, double[] v2) {
        return new double[] {
            v1[0] - v2[0],
            v1[1] - v2[1],
            v1[2] - v2[2]
        };
    }

    /**
     * Computes the torsion angle of 4 consecutive points.
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     * @return
     */
    public static double torsionAngle(double[] v1, double[] v2, double[] v3, double[] v4) {
        double[] ab = subtract(v1, v2);
        double[] cb = subtract(v3, v2);
        double[] bc = subtract(v2, v3);
        double[] dc = subtract(v4, v3);

        double[] abc = vectorProduct(ab, cb);
        double[] bcd = vectorProduct(bc, dc);

        double angle = angle(abc, bcd);
        /* calc the sign: */
        double[] vecprod = vectorProduct(abc, bcd);
        double val = dotProduct(cb, vecprod);
        if (val < 0.0) {
            angle = -angle;
        }

        return angle;
    }


    /**
     * Computes the vector product between 2 3D vectors.
     * @param v1 the first 3D vector
     * @param v2 the second 3D vector
     * @return the vector product
     */
    public static double[] vectorProduct(double[] v1, double[] v2) {
        return new double[] {
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        };
    }

    public static double distance14(double[] ca_p2, double[] ca_p1, double[] ca_tr, double[] ca_n1) {
        double d14 = distance(ca_p2, ca_n1);
        if(isLeftHandedTwist(ca_p2, ca_p1, ca_tr, ca_n1)) {
            d14 = -d14;
        }
        return d14;
    }

    private static boolean isLeftHandedTwist(double[] ca_p2, double[] ca_p1, double[] ca_tr, double[] ca_n1) {
        double[] d21 = subtract(ca_p1, ca_p2);
        double[] d31 = subtract(ca_tr, ca_p2);
        double[] d41 = subtract(ca_n1, ca_p2);

        return dotProduct(d21, vectorProduct(d31, d41)) < 0.0;
    }
}
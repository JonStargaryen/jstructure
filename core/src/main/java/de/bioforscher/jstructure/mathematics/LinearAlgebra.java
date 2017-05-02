package de.bioforscher.jstructure.mathematics;

/**
 * This class provides access to basic linear algebraic operations on either <code>double[]</code> or model instances
 * providing fields with e.g. coordinates.
 * Created by bittrich on 4/28/17.
 */
public class LinearAlgebra {
    private LinearAlgebra() {
        //TODO impl
    }

    public static PrimitiveDoubleArrayLinearAlgebra on(double[] value) {
        return new PrimitiveDoubleArrayLinearAlgebra(value);
    }

    /**
     * Caps a number to a defined interval.
     * @param lowerBound the lowest accepted value x_l
     * @param value the raw value
     * @param upperBound the highest accepted value x_u
     * @return the original value if it lies in the interval [x_l,x_u], else x_l or x_u
     */
    static double capToInterval(double lowerBound, double value, double upperBound) {
        return value < lowerBound ? lowerBound : value > upperBound ? upperBound : value;
    }

    public static class PrimitiveDoubleArrayLinearAlgebra {
        private double x, y, z;

        PrimitiveDoubleArrayLinearAlgebra(double[] value) {
            this.x = value[0];
            this.y = value[1];
            this.z = value[2];
        }

        public PrimitiveDoubleArrayLinearAlgebra add(double[] vectorToAdd) {
            this.x += vectorToAdd[0];
            this.y += vectorToAdd[1];
            this.z += vectorToAdd[2];
            return this;
        }

        public double angle(double[] secondVector) {
            double vDot = dotProduct(secondVector) / (norm() * LinearAlgebra.on(secondVector).norm());
            return Math.toDegrees(Math.acos(capToInterval(-1, vDot, 1)));
        }

        public double distance(double[] secondVector) {
            return Math.sqrt(distanceFast(secondVector));
        }

        public double distanceFast(double[] secondVector) {
            return (x - secondVector[0]) * (x - secondVector[0]) +
                    (y - secondVector[1]) * (y - secondVector[1]) +
                    (z - secondVector[2]) * (z - secondVector[2]);
        }

        public PrimitiveDoubleArrayLinearAlgebra divide(double divisor) {
            this.x /= divisor;
            this.y /= divisor;
            this.z /= divisor;
            return this;
        }

        public double dotProduct(double[] secondVector) {
            return x * secondVector[0] +
                    y * secondVector[1] +
                    z * secondVector[2];
        }

        public double[] getValue() {
            return new double[] { x, y, z };
        }

        public PrimitiveDoubleArrayLinearAlgebra multiply(double scalar) {
            this.x *= scalar;
            this.y *= scalar;
            this.z *= scalar;
            return this;
        }

        public PrimitiveDoubleArrayLinearAlgebra multiply(double[][] matrix) {
            double tx = matrix[0][0] * x + matrix[1][0] * y + matrix[2][0] * z;
            double ty = matrix[0][1] * x + matrix[1][1] * y + matrix[2][1] * z;
            double tz = matrix[0][2] * x + matrix[1][2] * y + matrix[2][2] * z;
            this.x = tx;
            this.y = ty;
            this.z = tz;
            return this;
        }

        public double norm() {
            return Math.sqrt(dotProduct(getValue()));
        }

        public PrimitiveDoubleArrayLinearAlgebra normalize() {
            divide(norm());
            return this;
        }

        public PrimitiveDoubleArrayLinearAlgebra substract(double[] secondVector) {
            this.x -= secondVector[0];
            this.y -= secondVector[1];
            this.z -= secondVector[2];
            return this;
        }

        public PrimitiveDoubleArrayLinearAlgebra vectorProduct(double[] secondVector) {
            double tx = y * secondVector[2] - z * secondVector[1];
            double ty = z * secondVector[0] - x * secondVector[2];
            double tz = x * secondVector[1] - y * secondVector[0];
            this.x = tx;
            this.y = ty;
            this.z = tz;
            return this;
        }
    }
}

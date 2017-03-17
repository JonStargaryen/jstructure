package de.bioforscher.explorer.membrane.o.dmodel;

/**
 * Conversions between raw continuous and discrete values and the RGB/HSV color space for front-end visualization.
 * Created by bittrich on 3/14/17.
 */
public class ColorMapping {
    public static double DEFAULT_S = 60.0;
    public static double DEFAULT_V = 70.0;
    public static double[] HSV_MISSING_VALUE = { 0.0, DEFAULT_S, DEFAULT_V };

    /**
     *
     * @param value in the interval [0.0, 1.0]
     * @return the hsv value
     */
    public static double[] continuousToHsv(double value) {
        return new double[] { value * 120.0, DEFAULT_S, DEFAULT_V };
    }

    public static double[] discreteToHsv(double value, int numberOfStates) {
        return continuousToHsv(value / (double) numberOfStates);
    }

    public static double[] discreteToHsv(int index, int numberOfStates) {
        return discreteToHsv((double) index, numberOfStates);
    }
}

package de.bioforscher.jstructure.model.structure;

/**
 * A model instance which provides atomic coordinates.
 * Created by bittrich on 4/28/17.
 */
public interface CoordinateProvider {
    /**
     * The spatial coordinates of this entity.
     * @return a <code>double[]</code> of dimensionality 3 resembling this entities spatial position
     */
    double[] getCoordinates();

    /**
     * Updates this entities spatial coordinates.
     * @param coordinates the new value, ought to be of dimensionality 3
     */
    void setCoordinates(double[] coordinates);
}

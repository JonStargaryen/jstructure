package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Atom;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Transformation extends FeatureContainerEntry {
    private static final Logger logger = LoggerFactory.getLogger(Transformation.class);
    private final double[] translation;
    private final double[][] rotation;

    /**
     * Rototranslation.
     * @param translation the translation vector
     * @param rotation the rotation matrix
     */
    public Transformation(double[] translation, double[][] rotation) {
        super(null);
        if(translation == null && rotation == null) {
            throw new IllegalArgumentException("cannot create transformation instance for missing translation and rotation");
        }
        this.translation = translation;
        this.rotation = rotation;
    }

    /**
     * Mere shift.
     * @param translation the translation vector
     */
    public Transformation(double[] translation) {
        this(translation, null);
    }

    /**
     * Mere rotation.
     * @param rotation the rotation matrix
     */
    public Transformation(double[][] rotation) {
        this(null, rotation);
    }

    public Atom transformCoordinates(Atom atom) {
        double[] vector = atom.getCoordinates();
        logger.trace("initial atom {}", atom.getPdbRepresentation());

        // apply rotation if needed
        if(rotation != null) {
            atom.setCoordinates(LinearAlgebra.on(vector).multiply(rotation).getValue());
        }

        vector = atom.getCoordinates();
        // apply transformation if needed
        if(translation != null) {
            atom.setCoordinates(LinearAlgebra.on(vector).add(translation).getValue());
        }

        logger.trace("transf. atom {}", atom.getPdbRepresentation());
        return atom;
    }

    public double[] getTranslation() {
        return translation;
    }

    public double[][] getRotation() {
        return rotation;
    }
}

package de.bioforscher.kinks;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;

import java.util.StringJoiner;

/**
 * Compute several features and compose arff.
 * Created by bittrich on 2/14/17.
 */
public class S02_ComposeKinkFinderArff {
    //TODO resolve in other packages
    private static final AbstractFeatureProvider asaCalculator = new AccessibleSurfaceAreaCalculator();
    private static final AbstractFeatureProvider plipAnnotator = new PLIPAnnotator();

    public static void main(String[] args) {
        StringJoiner header = new StringJoiner(System.lineSeparator());
        /*
         * @RELATION iris
         * @ATTRIBUTE sepallength  NUMERIC
         * @ATTRIBUTE sepalwidth   NUMERIC
         * @ATTRIBUTE petallength  NUMERIC
         * @ATTRIBUTE petalwidth   NUMERIC
         * @ATTRIBUTE class        {Iris-setosa,Iris-versicolor,Iris-virginica}
         */
        header.add("@RELATION kinks      NUMERIC");
        header.add("@ATTRIBUTE pdbId     NUMERIC");
        header.add("@ATTRIBUTE chainId   NUMERIC");
        header.add("@ATTRIBUTE start     NUMERIC");
        header.add("@ATTRIBUTE end       NUMERIC");
        header.add("@ATTRIBUTE kink      NUMERIC");
        header.add("@ATTRIBUTE gly       NUMERIC");
        header.add("@ATTRIBUTE pro       NUMERIC");
        header.add("@ATTRIBUTE rasaL4    NUMERIC");
        header.add("@ATTRIBUTE rasaL3    NUMERIC");
        header.add("@ATTRIBUTE rasaL2    NUMERIC");
        header.add("@ATTRIBUTE rasaL1    NUMERIC");
        header.add("@ATTRIBUTE rasa      NUMERIC");
        header.add("@ATTRIBUTE rasaR1    NUMERIC");
        header.add("@ATTRIBUTE rasaR2    NUMERIC");
        header.add("@ATTRIBUTE rasaR3    NUMERIC");
        header.add("@ATTRIBUTE rasaR4    NUMERIC");
        header.add("@ATTRIBUTE plipL4    NUMERIC");
        header.add("@ATTRIBUTE plipL3    NUMERIC");
        header.add("@ATTRIBUTE plipL2    NUMERIC");
        header.add("@ATTRIBUTE plipL1    NUMERIC");
        header.add("@ATTRIBUTE plip      NUMERIC");
        header.add("@ATTRIBUTE plipR1    NUMERIC");
        header.add("@ATTRIBUTE plipR2    NUMERIC");
        header.add("@ATTRIBUTE plipR3    NUMERIC");
        header.add("@ATTRIBUTE plipR4    NUMERIC");
        header.add("@ATTRIBUTE error     NUMERIC");
        header.add("@ATTRIBUTE angle     NUMERIC");

        /*
         * @DATA
         * 5.1,3.5,1.4,0.2,Iris-setosa
         * 4.9,3.0,1.4,0.2,Iris-setosa
         * 4.7,3.2,1.3,0.2,Iris-setosa
         * 4.6,3.1,1.5,0.2,Iris-setosa
         * 5.0,3.6,1.4,0.2,Iris-setosa
         * 5.4,3.9,1.7,0.4,Iris-setosa
         * 4.6,3.4,1.4,0.3,Iris-setosa
         * 5.0,3.4,1.5,0.2,Iris-setosa
         * 4.4,2.9,1.4,0.2,Iris-setosa
         * 4.9,3.1,1.5,0.1,Iris-setosa
         */
        header.add("@DATA");
    }
}

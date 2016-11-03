package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;

import java.util.ArrayList;
import java.util.List;

/**
 * The container for membranes.
 * Created by S on 03.11.2016.
 */

public class Membrane {
    final List<double[]> membraneAtoms;
    double[] point;
    double[] planePoint1;
    double[] planePoint2;
    double hphob;
    double hphil;
    double total;
    double qmax;
    public double[] normalVector;

    public Membrane(double[] c1, double[] c2, double hphob, double hphil) {
        this.planePoint1 = c1;
        this.planePoint2 = c2;
        this.hphob = hphob;
        this.hphil = hphil;
        this.total = hphob + hphil;
        this.normalVector = LinearAlgebra3D.subtract(c1, c2);
        this.membraneAtoms = new ArrayList<>();
    }

    public Membrane(double[] point, double[] c1, double[] c2, double hphob, double hphil, double total, double qmax) {
        this(c1, c2, hphob, hphil);
        this.point = point;
        this.total = total;
        this.qmax = qmax;
    }
}
package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.SelectionException;

/**
 * Represents torsion angles (phi, psi and omega) between 2 amino acids.
 * Created by bittrich on 5/23/17.
 */
public class TorsionAngles {
    private final double phi;
    private final double psi;
    private final double omega;

    public TorsionAngles(AminoAcid group1, AminoAcid group2) {
        try {
            Atom n1 = group1.getN();
            Atom ca1 = group1.getCa();
            Atom c1 = group1.getC();

            Atom n2 = group2.getN();
            Atom ca2 = group2.getCa();
            Atom c2 = group2.getC();

            this.phi = torsionAngle(c1, n2, ca2, c2);
            this.psi = torsionAngle(n1, ca1, c1, n2);
            this.omega = torsionAngle(ca1, c1, n2, ca2);
        } catch (NullPointerException e) {
            throw new SelectionException("missing backbone atoms for torsion angle calculation in " + group1 + " or " + group2);
        }
    }

    /**
     * Computes the torsion angle of 4 consecutive points.
     */
    private static double torsionAngle(Atom a1, Atom a2, Atom a3, Atom a4) {
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra ab = LinearAlgebra.on(a1).subtract(a2.getCoordinates());
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra cb = LinearAlgebra.on(a3).subtract(a2.getCoordinates());
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra bc = LinearAlgebra.on(a2).subtract(a3.getCoordinates());
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra dc = LinearAlgebra.on(a4).subtract(a3.getCoordinates());

        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra abc = ab.vectorProduct(cb);
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra bcd = bc.vectorProduct(dc);

        double angle = abc.angle(bcd);
        /* calc the sign: */
        double[] vecprod = abc.vectorProduct(bcd).getValue();
        double val = cb.dotProduct(vecprod);
        if (val < 0.0) {
            angle = -angle;
        }

        return angle;
    }

    public TorsionAngles(double phi, double psi, double omega) {
        this.phi = phi;
        this.psi = psi;
        this.omega = omega;
    }

    public double getPhi() {
        return phi;
    }

    public double getPsi() {
        return psi;
    }

    public double getOmega() {
        return omega;
    }
}

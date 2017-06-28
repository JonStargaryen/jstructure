package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;

import java.text.DecimalFormat;

/**
 * The element dedicated to describe secondary structure elements assigned by ASSP.
 * Created by bittrich on 6/28/17.
 */
public class ASSPSecondaryStructure extends GenericSecondaryStructure {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.##");
    private double t;
    private double h;
    private double vtor;
    private double r;
    private double[] u;
    private double ba;
    private boolean helicalParameters;
    // marks this as continuous stretch
    private boolean continuous;
    // twist, rise and radius fulfill alpha-helix or PPII (A, a, S, P, U)
    private String a1;
    // twist, rise and radius fulfill 3-10 (G, g, U)
    private String a2;
    // twist, rise and radius fulfill pi-helix (I, i, U)
    private String a3;
    private boolean helicalCharacteristics;

    ASSPSecondaryStructure(AbstractFeatureProvider featureProvider, SecondaryStructureElement secondaryStructure) {
        super(featureProvider, secondaryStructure);
    }

    void setHelicalParameters(double t, double h, double vtor, double r, double[] u) {
        this.t = t;
        this.h = h;
        this.vtor = vtor;
        this.r = r;
        this.u = u;
        this.helicalParameters = true;
    }

    void setCharacteristics(String a1, String a2, String a3) {
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;
        this.helicalCharacteristics = true;
    }

    void setBa(double ba) {
        this.ba = ba;
    }

    /**
     * The twist of a given fragment.
     * @return a double value in degree
     */
    public double getT() {
        return t;
    }

    /**
     * The height respectively rise of a fragment.
     * @return a double value in Angstrom
     */
    public double getH() {
        return h;
    }

    /**
     * The virtual torsion angle of a fragment.
     * @return a double value in degree
     */
    public double getVtor() {
        return vtor;
    }

    /**
     * The radius of a fragment.
     * @return a double value in Angstrom
     */
    public double getR() {
        return r;
    }

    double[] getU() {
        return u;
    }

    @Override
    public String toString() {
        if(helicalParameters) {
            if(continuous) {
                return "ASSP{" + (helicalCharacteristics ? a1 + " " + a2 + " " + a3 : "? ? ?") +
                        " Twist=" + DECIMAL_FORMAT.format(t) +
                        "°, h=" + DECIMAL_FORMAT.format(h) +
                        "A, Vtor=" + DECIMAL_FORMAT.format(vtor) +
                        "°, BA=" + DECIMAL_FORMAT.format(ba) +
                        "°, Radius=" + DECIMAL_FORMAT.format(r) +
                        "A}";
            } else {
                return "ASSP{non continuous}";
            }
        } else {
            return "ASSP{not assigned}";
        }
    }

    void setContinuous(boolean continuous) {
        this.continuous = continuous;
    }

    public boolean isContinuous() {
        return continuous;
    }
}

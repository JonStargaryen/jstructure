package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/**
 * The element dedicated to describe secondary structure elements assigned by ASSP.
 * Created by bittrich on 6/28/17.
 */
public class ASSPSecondaryStructure extends GenericSecondaryStructure {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.##", DecimalFormatSymbols.getInstance(Locale.US));
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
    private int stretchId;
    private String finalCharacteristic;

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

    void setCharacteristics(String[] helicalCharacteristics) {
        this.a1 = helicalCharacteristics[0];
        this.a2 = helicalCharacteristics[1];
        this.a3 = helicalCharacteristics[2];
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
                return "ASSP{" + (helicalCharacteristics ? a1 + " " + a2 + " " + a3 + " " + stretchId + "\t" : "? ? ? ?\t") +
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

    public int getStretchId() {
        return stretchId;
    }

    void setStretchId(int stretchId) {
        this.stretchId = stretchId;
    }

    boolean isUnassigned() {
        return "U".equals(a1) && "U".equals(a2) && "U".equals(a3);
    }

    String getAlpha() {
        return a1;
    }

    String getThree() {
        return a2;
    }

    String getPi() {
        return a3;
    }

    void setFinalCharacteristic(String finalCharacteristic) {
        this.finalCharacteristic = finalCharacteristic;
    }

    String getFinalCharacteristic() {
        return finalCharacteristic;
    }
}

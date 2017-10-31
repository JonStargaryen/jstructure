package de.bioforscher.jstructure.feature.sse.dssp;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureType;
import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

import java.util.Arrays;

/**
 * This class extends the basic container for secondary structure annotation,
 * including all the information used in the DSSP algorithm.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
@DefaultFeatureProvider(DictionaryOfProteinSecondaryStructure.class)
public class DSSPSecondaryStructure extends GenericSecondaryStructure {
    private double phi;
    private double psi;
    private double omega;
    private float kappa;

    private HBond accept1; // from CO of partner to NH of this
    private HBond accept2; // this is the donor of accept partner
    private HBond donor1; // from CO of this to NH of partner
    private HBond donor2; // this is the acceptor of donor partner

    // Symbols: starting '>', ending '<', or both 'X'.
    // Number means bracketed n-turn getResidue without h-bond
    private char[] turn;
    private boolean bend;

    private BetaBridge bridge1;
    private BetaBridge bridge2;

    DSSPSecondaryStructure(FeatureProvider featureProvider, SecondaryStructureType ss) {
        super(featureProvider, ss);

        this.phi = 360;
        this.psi = 360;
        this.omega = 360;

        this.accept1 = new HBond();
        this.accept2 = new HBond();
        this.donor1 = new HBond();
        this.donor2 = new HBond();

        this.bridge1 = null;
        this.bridge2 = null;

        this.turn = new char[] { ' ', ' ', ' ' };

        this.bend = false;
        this.kappa = 360;
    }

    @Override
    public String toString() {
        return "DSSPSecondaryStructure{" +
                "secondaryStructure=" + secondaryStructure +
                ", phi=" + phi +
                ", psi=" + psi +
                ", omega=" + omega +
                ", kappa=" + kappa +
                ", accept1=" + accept1 +
                ", accept2=" + accept2 +
                ", donor1=" + donor1 +
                ", donor2=" + donor2 +
                ", turn=" + Arrays.toString(turn) +
                ", bend=" + bend +
                ", bridge1=" + bridge1 +
                ", bridge2=" + bridge2 +
                '}';
    }

    public boolean isBend() {
        return this.bend;
    }

    public void setBend(boolean bend) {
        this.bend = bend;
    }

    public float getKappa() {
        return this.kappa;
    }

    public void setKappa(float kappa) {
        this.kappa = kappa;
    }

    public char[] getTurn() {
        return this.turn;
    }

    /**
     * Set the turn column corresponding to 3,4 or 5 helix patterns. If starting
     * > or ending < was set and the opposite is being set, the value will be
     * converted to X. If a number was set, it will be overwritten by the new
     * character.
     *
     * @param c
     *            character in the column
     * @param t
     *            turn of the helix {3,4,5}
     */
    public void setTurn(char c, int t) {
        if (this.turn[t - 3] == 'X')
            return;
        else if (this.turn[t - 3] == '<' && c == '>' || this.turn[t - 3] == '>' && c == '<') {
            this.turn[t - 3] = 'X';
        } else if (this.turn[t - 3] == '<' || this.turn[t - 3] == '>') {
            return;
        } else {
            this.turn[t - 3] = c;
        }
    }

    public HBond getAccept1() {
        return this.accept1;
    }

    public void setAccept1(HBond accept1) {
        this.accept1 = accept1;
    }

    public HBond getAccept2() {
        return this.accept2;
    }

    public void setAccept2(HBond accept2) {
        this.accept2 = accept2;
    }

    public HBond getDonor1() {
        return donor1;
    }

    public void setDonor1(HBond donor1) {
        this.donor1 = donor1;
    }

    public HBond getDonor2() {
        return this.donor2;
    }

    public void setDonor2(HBond donor2) {
        this.donor2 = donor2;
    }

    public double getPhi() {
        return this.phi;
    }

    public void setPhi(double phi) {
        this.phi = phi;
    }

    public double getPsi() {
        return this.psi;
    }

    public void setPsi(double psi) {
        this.psi = psi;
    }

    public double getOmega() {
        return this.omega;
    }

    public void setOmega(double omega) {
        this.omega = omega;
    }

    public BetaBridge getBridge1() {
        return this.bridge1;
    }

    public BetaBridge getBridge2() {
        return this.bridge2;
    }

    /**
     * Adds a Bridge to the getResidue. Each residue can only store two bridges. If
     * the residue contains already two Bridges, the Bridge will not be added
     * and the method returns false.
     *
     * @param bridge
     * @return false if the Bridge was not added, true otherwise
     */
    public boolean addBridge(BetaBridge bridge) {
        if (this.bridge1 == null) {
            this.bridge1 = bridge;
            return true;
        } else if (this.bridge1.equals(bridge)) {
            return true;
        } else if (this.bridge2 == null) {
            this.bridge2 = bridge;
            return true;
        } else if (this.bridge2.equals(bridge)) {
            return true;
        }
        // no space left, cannot add the bridge
        return false;
    }

    public void setBridge1(BetaBridge bridge1) {
        this.bridge1 = bridge1;
    }

    public void setBridge2(BetaBridge bridge2) {
        this.bridge2 = bridge2;
    }
}
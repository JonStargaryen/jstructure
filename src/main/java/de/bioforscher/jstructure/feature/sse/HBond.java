package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.model.structure.Residue;

/**
 * Container that represents a hydrogen bond. It contains the energy of the bond
 * in cal/mol and the partner index.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class HBond {
    private double energy;
    private Residue partner;

    public HBond() {
        this.energy = 0;
        this.partner = null;
    }

    public HBond(HBond o) {
        this.energy = o.energy;
        this.partner = o.partner;
    }

    @Override
    public HBond clone() {
        return new HBond(this);
    }

    @Override
    public String toString() {
        return this.partner + " | " + (this.energy / 1000.0);
    }

    public double getEnergy() {
        return this.energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public Residue getPartner() {
        return this.partner;
    }

    public void setPartner(Residue partner) {
        this.partner = partner;
    }
}
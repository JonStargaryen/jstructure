package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.prototype.AbstractGroup;

/**
 * The abstract representation of an amino acid.
 * Created by bittrich on 5/24/17.
 */
public abstract class AminoAcid extends AbstractGroup implements StandardAminoAcidIndicator {
    private Atom n;
    private Atom ca;
    private Atom c;
    private Atom o;
    private Atom h;

    public Atom getN() {
        return n;
    }

    public Atom getCa() {
        return ca;
    }

    public Atom getC() {
        return c;
    }

    public Atom getO() {
        return o;
    }

    public Atom getH() {
        return h;
    }
}

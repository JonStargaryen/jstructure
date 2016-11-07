package de.bioforscher.jstructure.model.structure;

import java.util.Optional;

/**
 *
 * Created by S on 06.10.2016.
 */
public class Nucleotide extends Group {
    public Nucleotide(String pdbName, int residueNumber) {
        super(pdbName, residueNumber);
    }
//TODO implement
    public Optional<Atom> O3Prime() {
        return null;
    }

    public Optional<Atom> O5Prime() {
        return null;
    }

    public Optional<Atom> Phosphate() {
        return null;
    }
}

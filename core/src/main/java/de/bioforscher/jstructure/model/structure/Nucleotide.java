package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;

import java.util.Optional;

/**
 * Created by S on 08.11.2016.
 */
public class Nucleotide extends Group {
    public Nucleotide(String pdbName, int residueNumber) {
        super(pdbName, residueNumber);
    }

    public Atom getO3Prime() {
        return getAtomByName(AtomNameFilter.O3PRIME_FILTER);
    }

    public Optional<Atom> findO3Prime() {
        return tryToGetAtomByName(AtomNameFilter.O3PRIME_FILTER);
    }

    public Atom getO5Prime() {
        return getAtomByName(AtomNameFilter.O5PRIME_FILTER);
    }

    public Optional<Atom> findO5Prime() {
        return tryToGetAtomByName(AtomNameFilter.O5PRIME_FILTER);
    }

    public Atom getPhosphate() {
        return getAtomByName(AtomNameFilter.PHOSPHATE_FILTER);
    }

    public Optional<Atom> findPhosphate() {
        return tryToGetAtomByName(AtomNameFilter.PHOSPHATE_FILTER);
    }
}

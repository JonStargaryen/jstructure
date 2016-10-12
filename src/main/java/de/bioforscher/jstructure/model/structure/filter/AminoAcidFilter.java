package de.bioforscher.jstructure.model.structure.filter;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Residue;

import java.util.List;
import java.util.Objects;
import java.util.function.Predicate;

/**
 * Tests whether an atom is of a certain <tt>PDB</tt> name. This allows for filtering for CA or backbone atoms.
 */
public class AminoAcidFilter implements Predicate<Residue> {
    private final List<AminoAcid> acceptedAminoAcids;

    /**
     * Construct an amino acid specific filter.
     * @param acceptedAminoAcids references to the amino acids which will be accepted by this filter
     */
    public AminoAcidFilter(List<AminoAcid> acceptedAminoAcids) {
        this.acceptedAminoAcids = Objects.requireNonNull(acceptedAminoAcids);
    }

    @Override
    public boolean test(Residue residue) {
        return this.acceptedAminoAcids.contains(residue.getAminoAcid());
    }
}

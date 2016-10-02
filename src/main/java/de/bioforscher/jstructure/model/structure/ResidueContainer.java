package de.bioforscher.jstructure.model.structure;

import java.util.List;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a residue container (mostly a {@link Chain}).
 * Created by S on 30.09.2016.
 */
public interface ResidueContainer extends AtomContainer {
    /**
     * Access to all residues associated to this container.
     * @return a stream of residues
     */
    Stream<Residue> residues();

    /**
     * Tests whether an atom is of a certain <tt>PDB</tt> name. This allows for filtering for CA or backbone atoms.
     */
    class AminoAcidFilter implements Predicate<Residue> {
        private List<AminoAcid> acceptedAminoAcids;
        private boolean negate;

        /**
         * Construct an amino acid specific filter.
         * @param acceptedAminoAcids references to the amino acids which will be accepted by this filter
         * @param negate when true the given amino acids will evaluate to <code>false</code> and everything else will
         *               evaluate to <code>true</code>
         * @see Stream#filter(Predicate)
         */
        public AminoAcidFilter(List<AminoAcid> acceptedAminoAcids, boolean negate) {
            Objects.requireNonNull(acceptedAminoAcids);
            this.acceptedAminoAcids = acceptedAminoAcids;
            this.negate = negate;
        }

        @Override
        public boolean test(Residue residue) {
            return negate != this.acceptedAminoAcids.contains(residue.getAminoAcid());
        }

        /**
         * Creates a negated instance of this filter. Negating twice will result in a filter with the initial behaviour.
         * @return a negated instance of this filter
         */
        public AminoAcidFilter getNegatedInstance() {
            return new AminoAcidFilter(acceptedAminoAcids, !negate);
        }
    }
}
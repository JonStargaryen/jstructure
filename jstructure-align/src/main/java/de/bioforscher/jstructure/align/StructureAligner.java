package de.bioforscher.jstructure.align;

/**
 * Specification of structure alignment algorithms. Queries are build via {@link StructureAlignmentBuilder}.
 * Created by bittrich on 6/19/17.
 */
public interface StructureAligner {
    /**
     * Structurally align 2 containers of atoms.
     * @param builder the builder providing all required settings
     * @return a structure alignment by the given criteria
     * @throws AlignmentException when no alignment could be achieved
     */
    StructureAlignmentResult align(StructureAlignmentBuilder.StructureAlignmentStep builder) throws AlignmentException;
}

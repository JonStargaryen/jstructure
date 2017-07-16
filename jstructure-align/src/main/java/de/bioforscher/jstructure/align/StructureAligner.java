package de.bioforscher.jstructure.align;

/**
 * Specification of structure alignment algorithms. Queries are build via {@link StructureAlignmentQuery}.
 * Created by bittrich on 6/19/17.
 */
public interface StructureAligner {
    /**
     * Structurally align 2 containers of atoms.
     * @param structureAlignmentQuery the query providing all required settings
     * @return a structure alignment by the given criteria
     * @throws AlignmentException when no alignment could be achieved
     */
    StructureAlignmentResult align(StructureAlignmentQuery structureAlignmentQuery) throws AlignmentException;
}

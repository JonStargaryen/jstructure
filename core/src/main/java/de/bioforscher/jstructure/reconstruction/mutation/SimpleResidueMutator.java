package de.bioforscher.jstructure.reconstruction.mutation;

import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.alignment.StructureAlignmentResult;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.parser.CIFParser;

/**
 * Naive approach to mutate a residue: look for the maximum agreement between original and new amino acid and superimpose
 * new one onto overlapping region.
 * Created by bittrich on 4/21/17.
 */
public class SimpleResidueMutator {
    public void mutatePosition(Group group, AminoAcidFamily targetAminoAcid) {
        //TODO this needs the amino-acid-scaffold implementation
        Group scaffold = (Group) group.getParentChain()
                .select()
                .aminoAcids(targetAminoAcid)
                .asGroup()
                .createCopy();

        StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(group, scaffold);
        alignmentResult.transform(scaffold);

        // rename amino acid as this determines the newly place entity
        group.setGroupInformation(CIFParser.parseLigandInformation(targetAminoAcid.getThreeLetterCode()));

        group.clearAtoms();
        scaffold.atoms().forEach(group::addAtom);

        //TODO abstract interface, renumber atoms
    }
}

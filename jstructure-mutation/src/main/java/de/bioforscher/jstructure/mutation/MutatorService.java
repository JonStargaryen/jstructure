package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

/**
 * Introduces a mutation for a given residue in a protein structure.
 *
 * Assumptions made:
 * <ul>
 *     <li>all chains are considered individually and isolated from potentially other chains present in a protein</li>
 * </ul>
 * Created by bittrich on 7/12/17.
 */
public interface MutatorService {
    /**
     * Mutate a given position in a protein.
     * @param originalProtein the original protein which will not be manipulated
     * @param aminoAcidToMutate the amino acid to mutate
     * @param targetAminoAcid the desired target amino acid
     * @return a new {@link Protein} instance featuring the mutation
     */
    Protein mutateAminoAcid(Protein originalProtein, AminoAcid aminoAcidToMutate, AminoAcid.Family targetAminoAcid);
}

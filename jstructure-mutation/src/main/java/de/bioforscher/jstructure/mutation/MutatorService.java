package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Structure;
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
     * @param chainIdentifier the chain where the mutation shall be introduced
     * @param residueToMutate the amino acid to mutate
     * @param targetAminoAcid the desired target amino acid
     * @return a new {@link Structure} instance featuring the mutation
     */
    Structure mutateAminoAcid(Structure originalProtein,
                              ChainIdentifier chainIdentifier,
                              ResidueIdentifier residueToMutate,
                              AminoAcid.Family targetAminoAcid);
}

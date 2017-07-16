package de.bioforscher.jstructure.mutation.old;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

/**
 * Introduces a mutation for a given residue in a protein structure.
 * Created by bittrich on 7/12/17.
 */
@Deprecated
public interface ResidueMutatorService {
    /**
     * Mutate a given position in a protein.
     * @param originalProtein the original protein which will not be manipulated
     * @param chainId the chainId of the mutation
     * @param residueNumber the residueNumber of the mutation
     * @param targetAminoAcid the desired target amino acid
     * @return a new {@link Structure} instance featuring the mutation
     */
    Structure mutateResidue(Structure originalProtein, String chainId, int residueNumber, AminoAcid.Family targetAminoAcid);
}

package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

/**
 * Screens for ligands around a given amino acid.
 * Created by bittrich on 7/13/17.
 */
public interface LigandContactScreener {
    /**
     * Determines the number of ligands in spatial proximity of an amino acid.
     * @param protein the container to search in
     * @param aminoAcid the peculiar amino acid
     * @return the number of observed contacts within a given cutoff
     */
    int determineNumberOfLigandContacts(Protein protein, AminoAcid aminoAcid);
}

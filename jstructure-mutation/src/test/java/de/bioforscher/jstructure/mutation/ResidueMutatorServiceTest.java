package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Test;

/**
 * Test for the residue mutator.
 * Created by bittrich on 7/6/17.
 */
public class ResidueMutatorServiceTest {
    @Test
    public void shouldIntroduceMutation() {
        Protein protein = ProteinParser.source("2lzm").parse();
        Protein mutatedProtein = new ResidueMutatorService().mutate(protein, "A", 27, AminoAcid.Family.ALANINE);

        System.out.println(protein.getPdbRepresentation());
        System.out.println(mutatedProtein.getPdbRepresentation());
    }
}
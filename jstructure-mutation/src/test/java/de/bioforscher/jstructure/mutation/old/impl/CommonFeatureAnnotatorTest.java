package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.old.ResidueMutatorService;
import org.junit.Before;
import org.junit.Test;

/**
 * Test integration of {@link ResidueMutatorService} and {@link CommonFeatureAnnotator}.
 * Created by bittrich on 7/13/17.
 */
public class CommonFeatureAnnotatorTest {
    private ResidueMutatorServiceImpl residueMutatorService;

    @Before
    public void setup() {
        residueMutatorService = new ResidueMutatorServiceImpl();
    }

    @Test
    public void shouldAnnotateMutatedProtein() {
        String chainId = "A";
        int position = 73;
        Structure protein = StructureParser.source("1aqe")
                .minimalParsing(true)
                .parse();
        Structure mutatedProtein = residueMutatorService.mutateResidue(protein, chainId, position, AminoAcid.Family.GLUTAMIC_ACID);
        CommonFeatureAnnotator.annotateProtein(mutatedProtein);

        AminoAcid mutatedAminoAcid = mutatedProtein.select()
                .chainName(chainId)
                .residueNumber(position)
                .asAminoAcid();
        System.out.println(new FeatureVector(mutatedAminoAcid));
    }
}
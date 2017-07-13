package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.ResidueMutatorService;
import org.junit.Before;
import org.junit.Test;

import java.nio.file.Paths;

/**
 * Test integration of {@link ResidueMutatorService} and {@link CommonFeatureAnnotator}.
 * Created by bittrich on 7/13/17.
 */
public class CommonFeatureAnnotatorTest {
    private ResidueMutatorServiceImpl residueMutatorService;

    @Before
    public void setup() {
        residueMutatorService = new ResidueMutatorServiceImpl();
        ProteinParser.OptionalSteps.setLocalPdbDirectory(Paths.get("/home/bittrich/pdb/"));
    }

    @Test
    public void shouldAnnotateMutatedProtein() {
        String chainId = "A";
        int position = 73;
        Protein protein = ProteinParser.localPdb("1aqe")
                .minimalParsing(true)
                .parse();
        Protein mutatedProtein = residueMutatorService.mutateResidue(protein, chainId, position, AminoAcid.Family.GLUTAMIC_ACID);
        CommonFeatureAnnotator.annotateProtein(mutatedProtein);

        AminoAcid mutatedAminoAcid = mutatedProtein.select()
                .chainName(chainId)
                .residueNumber(position)
                .asAminoAcid();
        System.out.println(new FeatureVector(mutatedAminoAcid));
    }
}
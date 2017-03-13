package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test the mapping between pdbId and UniProt.
 * Created by bittrich on 3/13/17.
 */
public class UniProtMapperTest {
    private AbstractFeatureProvider uniProtMapper = new UniProtMapper();

    @Test
    public void shouldAnnotate4LG6() {
        Protein protein = ProteinParser.source("4LG6").parse();
        uniProtMapper.process(protein);
        String uniProtIdA = protein.getChains().get(0).getFeature(String.class, SiftsParser.UNIPROT_ID);
        String uniProtIdB = protein.getChains().get(1).getFeature(String.class, SiftsParser.UNIPROT_ID);
        Assert.assertEquals(uniProtIdA, "Q9H9E1");
        Assert.assertEquals(uniProtIdB, "Q9H0W5");
    }
}

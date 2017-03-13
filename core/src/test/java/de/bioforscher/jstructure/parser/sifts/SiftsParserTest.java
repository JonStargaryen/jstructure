package de.bioforscher.jstructure.parser.sifts;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the {@link SiftsParser}.
 * Created by bittrich on 3/13/17.
 */
public class SiftsParserTest {
    private AbstractFeatureProvider featureProvider;

    @Before
    public void setup() {
        this.featureProvider = FeatureProviderRegistry.resolve(SiftsParser.EC_NUMBER);
    }

    @Test
    public void shouldMap102l() {
        String id = "102l";
        Protein protein = ProteinParser.source(id).parse();
        featureProvider.process(protein);

        Chain chain = protein.getChains().get(0);
        Assert.assertEquals("3.2.1.17", chain.getFeature(String.class, SiftsParser.EC_NUMBER));
        Assert.assertEquals("P00720", chain.getFeature(String.class, SiftsParser.UNIPROT_ID));
        Assert.assertEquals("PF00959", chain.getFeature(String.class, SiftsParser.PFAM_ID));
    }
}
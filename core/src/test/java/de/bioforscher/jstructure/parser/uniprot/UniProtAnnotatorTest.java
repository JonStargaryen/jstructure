package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test the UniProt annotator.
 * Created by bittrich on 3/2/17.
 */
public class UniProtAnnotatorTest {
    private AbstractFeatureProvider uniProtAnnotator = FeatureProviderRegistry.resolve(UniProtAnnotator.UNIPROT_ANNOTATION);
    @Test
    public void shouldAnnotate4LG6() {
        Protein protein = ProteinParser.parseProteinById("4LG6");
        uniProtAnnotator.process(protein);

        protein.chains().forEach(chain -> {
            String uniProtId = chain.getFeature(String.class, UniProtAnnotator.UNIPROT_ID);
            Assert.assertEquals(uniProtId, "Q9H0W5");
            UniProtAnnotationContainer container = chain.getFeature(UniProtAnnotationContainer.class, UniProtAnnotator.UNIPROT_ANNOTATION);
            Assert.assertTrue(container.getReferences().size() == 6);
            Assert.assertTrue(container.getMutagenesisSites().size() == 3);
            Assert.assertTrue(container.getNaturalVariants().size() == 0);
        });
    }
}
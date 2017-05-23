package de.bioforscher.jstructure.feature.uniprot;

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
    private AbstractFeatureProvider uniProtAnnotator = FeatureProviderRegistry.resolve(UniProtAnnotationContainer.class);

    @Test
    public void shouldAnnotate4LG6() {
        Protein protein = ProteinParser.source("4LG6").parse();
        uniProtAnnotator.process(protein);

        protein.chains().forEach(chain -> {
            UniProtAnnotationContainer container = chain.getFeatureContainer().getFeature(UniProtAnnotationContainer.class);
            Assert.assertTrue(container.getReferences().size() == 6);
            Assert.assertTrue(container.getMutagenesisSites().size() == 3);
            Assert.assertTrue(container.getNaturalVariants().size() == 0);
        });
    }
}
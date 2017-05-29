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
    public void shouldAnnotate4lgd() {
        Protein protein = ProteinParser.source("4lgd").parse();
        uniProtAnnotator.process(protein);

        protein.select()
                .chainName("A", "B", "C", "D")
                .asFilteredChains()
                .forEach(chain -> {
                        UniProtAnnotationContainer container = chain.getFeatureContainer().getFeature(UniProtAnnotationContainer.class);
        });
    }

    @Test
    public void shouldAnnotate4lg6() {
        Protein protein = ProteinParser.source("4lg6").parse();
        uniProtAnnotator.process(protein);

        protein.select()
                .chainName("A")
                .asFilteredChains()
                .forEach(chain -> {
                    UniProtAnnotationContainer container = chain.getFeatureContainer().getFeature(UniProtAnnotationContainer.class);
                    Assert.assertTrue(container.getReferences().size() == 6);
                    Assert.assertTrue(container.getMutagenesisSites().size() == 3);
                    Assert.assertTrue(container.getNaturalVariants().size() == 0);
                });
    }
}
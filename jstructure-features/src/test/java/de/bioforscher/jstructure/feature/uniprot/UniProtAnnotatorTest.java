package de.bioforscher.jstructure.feature.uniprot;

import de.bioforscher.jstructure.feature.uniprot.old.UniProtAnnotationContainer;
import de.bioforscher.jstructure.feature.uniprot.old.UniProtAnnotator;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test the UniProt annotator.
 * Created by bittrich on 3/2/17.
 */
@Deprecated
public class UniProtAnnotatorTest {
    private UniProtAnnotator uniProtAnnotator = new UniProtAnnotator();

    @Test
    public void shouldAnnotateNaturalVariants() {
        Protein protein = ProteinParser.source("4qei").parse();
        uniProtAnnotator.process(protein);
        UniProtAnnotationContainer container = protein.select()
                                                      .chainName("A")
                                                      .asChain()
                                                      .getFeatureContainer()
                                                      .getFeature(UniProtAnnotationContainer.class);
        Assert.assertTrue(container.getNaturalVariants().size() > 0);
    }

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
                    Assert.assertTrue(6 <= container.getReferences().size());
                    Assert.assertTrue(3 <= container.getMutagenesisSites().size());
                    Assert.assertTrue(0 <= container.getNaturalVariants().size());
                });
    }
}
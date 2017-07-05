package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

import java.io.IOException;

/**
 * Test the PLIP annotator.
 * Created by bittrich on 2/9/17.
 */
public class PLIPAnnotatorTest {
    private PLIPAnnotator plipAnnotator = new PLIPAnnotator();

    @Test
    public void shouldAnnotateSingleProtein() throws IOException {
        Protein protein = ProteinParser.source("4BPM").parse();
        plipAnnotator.process(protein);
        PLIPInteractionContainer plipInteractionContainer = protein.getFeatureContainer().getFeature(PLIPInteractionContainer.class);
        plipInteractionContainer.getInteractions().forEach(System.out::println);
    }

    @Test
    public void shouldAnnotateProteinWithWaterBridges() {
        Protein protein = ProteinParser.source("5A1S").parse();
        plipAnnotator.process(protein);
    }

    @Test
    public void shouldProteinWithIllReference() {
        Protein protein = ProteinParser.source("3AOU").parse();
        plipAnnotator.process(protein);
    }

    @Test
    public void shouldAnnotateProteinWithHalogenBond() {
        Protein protein = ProteinParser.source("4BPM").parse();
        plipAnnotator.process(protein);
    }
}
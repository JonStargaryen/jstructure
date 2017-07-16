package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
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
        Structure protein = StructureParser.source("4BPM").parse();
        plipAnnotator.process(protein);
        PLIPInteractionContainer plipInteractionContainer = protein.getFeatureContainer().getFeature(PLIPInteractionContainer.class);
        plipInteractionContainer.getInteractions().forEach(System.out::println);
    }

    @Test
    public void shouldAnnotateProteinWithWaterBridges() {
        Structure protein = StructureParser.source("5A1S").parse();
        plipAnnotator.process(protein);
    }

    @Test
    public void shouldProteinWithIllReference() {
        Structure protein = StructureParser.source("3AOU").parse();
        plipAnnotator.process(protein);
    }

    @Test
    public void shouldAnnotateProteinWithHalogenBond() {
        Structure protein = StructureParser.source("4BPM").parse();
        plipAnnotator.process(protein);
    }
}
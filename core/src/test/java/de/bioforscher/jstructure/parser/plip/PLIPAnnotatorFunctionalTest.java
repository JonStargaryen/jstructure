package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * Test the PLIP annotator.
 * Created by bittrich on 2/9/17.
 */
public class PLIPAnnotatorFunctionalTest {
    private AbstractFeatureProvider plipAnnotator;

    @Before
    public void setup() {
        plipAnnotator = new PLIPAnnotator();
    }

    @Test
    public void shouldAnnotateSingleProtein() throws IOException {
        Protein protein = ProteinParser.source("4BPM").parse();
        plipAnnotator.process(protein);
        PLIPInteractionContainer plipInteractionContainer = protein.getFeature(PLIPInteractionContainer.class, PLIPAnnotator.PLIP_INTERACTIONS);
        plipInteractionContainer.getInteractions().forEach(System.out::println);
    }

    @Test
    public void shouldAnnotateDataSet() throws IOException {
        System.err.println("skipping PLIP computation for data set");
//        final String listPath = "/home/bittrich/git/phd_sb_repo/data/pdbtm_alpha_nr.list.txt";
//        Files.lines(Paths.get(listPath))
//                .map(line -> line.substring(0, 4))
//                .map(id -> ProteinParser.source(id).parse())
//                .forEach(plipAnnotator::process);
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
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
        Protein protein = ProteinParser.parseProteinById("4BPM");
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
//                .map(ProteinParser::parseProteinById)
//                .forEach(plipAnnotator::process);
    }

    @Test
    public void shouldAnnotateProteinWithWaterBridges() {
        Protein protein = ProteinParser.parseProteinById("5A1S");
        plipAnnotator.process(protein);
    }

    @Test
    public void shouldProteinWithIllReference() {
        Protein protein = ProteinParser.parseProteinById("3AOU");
        plipAnnotator.process(protein);
    }

    @Test
    public void shouldAnnotateProteinWithHalogenBond() {
        Protein protein = ProteinParser.parseProteinById("4BPM");
        plipAnnotator.process(protein);
    }
}
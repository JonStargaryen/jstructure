package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Test the PLIP annotator.
 * Created by bittrich on 2/9/17.
 */
public class PLIPAnnotatorFunctionalTest {
    private AbstractFeatureProvider plipAnnotator;

    @Before
    public void setup() {
        plipAnnotator = FeatureProviderRegistry.resolve(PLIPAnnotator.PLIP_INTERACTIONS);
    }
    @Test
    public void shouldAnnotateDataSet() throws IOException {
        final String listPath = "/home/bittrich/git/phd_sb_repo/data/dataset/nrpdbtm/pdbtm_alpha_nr.list.txt";
        Files.lines(Paths.get(listPath))
                .map(line -> line.substring(0, 4))
                .map(ProteinParser::parseProteinById)
                .forEach(plipAnnotator::process);
    }
}
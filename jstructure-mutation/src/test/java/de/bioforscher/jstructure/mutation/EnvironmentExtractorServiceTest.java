package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests for the environment extractor service.
 * Created by bittrich on 7/6/17.
 */
public class EnvironmentExtractorServiceTest {
    private EnvironmentExtractorService environmentExtractor;
    private PLIPAnnotator plipAnnotator;

    @Before
    public void setup() {
        environmentExtractor = new EnvironmentExtractorService();
        plipAnnotator = new PLIPAnnotator();
    }

    @Test
    public void shouldExtractEnvironment() {
        Protein protein = ProteinParser.source("2lzm").parse();
        plipAnnotator.process(protein);
        Environment environment = environmentExtractor.extractEnvironment(protein, "A", 27);
        System.out.println(environment);
        environmentExtractor.searchForSimilarEnvironmentInDatabase(environment).forEach(System.out::println);
    }
}
package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.nodes.Document;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;

/**
 * Ensure the functionality of the {@link OrientationsOfProteinsInMembranesAnnotator}.
 * Created by bittrich on 5/17/17.
 */
public class OrientationsOfProteinsInMembranesAnnotatorTest {
    private OrientationsOfProteinsInMembranesAnnotator annotator = new OrientationsOfProteinsInMembranesAnnotator();

    @Test
    public void shouldMoveToRepresentative() throws IOException {
        Document document = annotator.getDocument(OrientationsOfProteinsInMembranesAnnotator.SEARCH_URL + "1brr");
        Assert.assertTrue(document.text().startsWith("1m0l"));
    }

    @Test
    public void shouldFetchResultsDirectly() throws IOException {
        Document document = annotator.getDocument(OrientationsOfProteinsInMembranesAnnotator.SEARCH_URL + "1m0l");
        Assert.assertTrue(document.text().startsWith("1m0l"));
    }

    @Test
    public void shouldHandleTransMembraneProtein() {

    }

    @Test(expected = ComputationException.class)
    public void shouldFailOnNonTransMembraneProtein() {
        annotator.process(ProteinParser.source("1pmm").parse());
    }
}
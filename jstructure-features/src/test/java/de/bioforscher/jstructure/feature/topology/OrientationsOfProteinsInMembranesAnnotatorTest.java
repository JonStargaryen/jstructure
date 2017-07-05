package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.nodes.Document;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * Ensure the functionality of the {@link OrientationsOfProteinsInMembranesAnnotator}.
 * Created by bittrich on 5/17/17.
 */
public class OrientationsOfProteinsInMembranesAnnotatorTest {
    private OrientationsOfProteinsInMembranesAnnotator annotator = new OrientationsOfProteinsInMembranesAnnotator();
    private Protein protein;

    @Before
    public void setup() {
        protein = ProteinParser.source("1brr").parse();
    }

    @Test(expected = ComputationException.class)
    public void shouldFailForNotPresentProtein() {
        OrientationsOfProteinsInMembranesAnnotator.getDocument("2ac6");
    }

    @Test
    public void shouldHandle2k21() {
        String id = "2k21";
        Protein protein = ProteinParser.source(id).parse();
        annotator.process(protein);
        MembraneContainer container = protein.getFeatureContainer().getFeature(MembraneContainer.class);
        Assert.assertEquals(0, container.getTransMembraneHelices().size());
        Assert.assertTrue(container.getMembraneAtoms().size() > 0);
    }

    @Test
    public void shouldHandle3abk() {
        String id = "3abk";
        Protein protein = ProteinParser.source(id).parse();
        Document document = OrientationsOfProteinsInMembranesAnnotator.getDocument(id);
        annotator.process(protein, document);
    }

    @Test
    public void shouldHandle4rz0() {
        annotator.process(ProteinParser.source("4zr0").parse());
    }

    @Test
    public void shouldMoveToRepresentative() throws IOException {
        Document document = OrientationsOfProteinsInMembranesAnnotator.getDocumentInternal(OrientationsOfProteinsInMembranesAnnotator.SEARCH_URL + "1brr");
        Assert.assertTrue(document.text().startsWith("1m0l"));
    }

    @Test
    public void shouldFetchResultsDirectly() throws IOException {
        Document document = OrientationsOfProteinsInMembranesAnnotator.getDocumentInternal(OrientationsOfProteinsInMembranesAnnotator.SEARCH_URL + "1m0l");
        Assert.assertTrue(document.text().startsWith("1m0l"));
    }

    @Test
    public void shouldHandleTransMembraneProtein() {
        annotator.process(protein);
    }

    @Test
    public void shouldHandleMalformedData() {
        // id: 5a1s - chain B misses bracket
        // B - Tilt: 10° - Segments: 1( 29- 44), 2( 51- 73), 3( 82- 94), 4( 116- 132), 5( 137- 166), 6 208- 230), 7( 265- 286), 8( 300- 314), 9( 322- 347), 10( 359- 384), 11( 422- 442)
        annotator.process(ProteinParser.source("5a1s").parse());
    }

    @Test(expected = ComputationException.class)
    public void shouldFailOnNonTransMembraneProtein() {
        annotator.process(ProteinParser.source("1pmm").parse());
    }
}
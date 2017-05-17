package de.bioforscher.jstructure.feature.topology;

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
        String id = "1brr";
        annotator.getDocument(OrientationsOfProteinsInMembranesAnnotator.SEARCH_URL + "1brr");
        //TODO impl
//        Assert.assertTrue("protein name should be different - found " + protein.getName() + ", expected not " + id, !protein.getName().equals(id));
    }

    @Test
    public void shouldFetchResultsDirectly() throws IOException {
        String id = "1m0l";
        annotator.getDocument(OrientationsOfProteinsInMembranesAnnotator.SEARCH_URL + "1m0l");
        //TODO impl
//        Assert.assertTrue("protein name should be " + id + " - found " + protein.getName(), protein.getName().equals(id));
    }
}
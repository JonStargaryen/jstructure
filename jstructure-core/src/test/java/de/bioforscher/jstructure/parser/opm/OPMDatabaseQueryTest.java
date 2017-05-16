package de.bioforscher.jstructure.parser.opm;

import de.bioforscher.jstructure.model.structure.Protein;
import org.junit.Assert;
import org.junit.Test;

/**
 * Ensure the functionality of the {@link OPMDatabaseQuery}.
 * Created by bittrich on 3/13/17.
 */
public class OPMDatabaseQueryTest {
    @Test
    public void shouldMoveToRepresentative() {
        String id = "1brr";
        Protein protein = OPMDatabaseQuery.parseAnnotatedProteinById(id);
        Assert.assertTrue("protein name should be different - found " + protein.getName() + ", expected not " + id, !protein.getName().equals(id));
    }

    @Test
    public void shouldFetchResultsDirectly() {
        String id = "1m0l";
        Protein protein = OPMDatabaseQuery.parseAnnotatedProteinById(id);
        Assert.assertTrue("protein name should be " + id + " - found " + protein.getName(), protein.getName().equals(id));
    }
}
package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

/**
 * Test for the ConSurf integration.
 * Created by bittrich on 7/19/17.
 */
public class ConSurfDBQueryTest {
    private ConSurfDBQuery conSurfDBQuery;

    @Before
    public void setup() {
        conSurfDBQuery = new ConSurfDBQuery();
    }

    @Test
    public void shouldQueryConSurf() {
        String pdbId = "1acj";
        Structure protein = StructureParser.source(pdbId)
                .minimalParsing(true)
                .parse();
        conSurfDBQuery.process(protein);

        protein.aminoAcids()
                .map(aminoAcid -> aminoAcid.getIdentifier() + " " + aminoAcid.getFeature(ConSurfDBQuery.ConSurfScore.class).getScore())
                .forEach(System.out::println);
    }

    @Test
    public void shouldParseResultFile() {
        conSurfDBQuery.parseGradesFiles(TestUtils.getResourceAsInputStream("consurf.grades"));
    }
}
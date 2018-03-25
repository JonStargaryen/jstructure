package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.query.StructureAlignmentQuery;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

import static de.bioforscher.testutil.TestUtils.SupportedProtein.PDB_1ACJ;
import static de.bioforscher.testutil.TestUtils.SupportedProtein.PDB_1BRR;

public class TMAlignServiceTest {
    private Structure reference;
    private Structure query;
    private TMAlignService instance;

    @Before
    public void setup() {
        reference = StructureParser.fromInputStream(TestUtils.getProteinInputStream(PDB_1BRR))
                .minimalParsing(true)
                .parse();
        query = StructureParser.fromInputStream(TestUtils.getProteinInputStream(PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        instance = TMAlignService.getInstance();
    }

    @Test
    public void shouldAlignStructures() {
        System.out.println(instance.process(StructureAlignmentQuery.builder()
                .reference(reference)
                .query(query)
                .build()));
    }
}
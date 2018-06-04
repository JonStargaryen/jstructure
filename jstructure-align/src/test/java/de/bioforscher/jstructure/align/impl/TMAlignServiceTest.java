package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.align.query.StructureAlignmentQuery;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
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
    @Ignore
    public void shouldAlignStructures() throws AlignmentException {
        TMAlignAlignmentResult tmAlignAlignmentResult = instance.process(StructureAlignmentQuery.builder()
                .reference(reference)
                .query(query)
                .build());

        Assert.assertEquals("rmsd does not match",
                6.36,
                tmAlignAlignmentResult.getRootMeanSquareDeviation().getScore(),
                TestUtils.TOLERANT_ERROR_MARGIN);
        Assert.assertEquals("tmscore does not match",
                0.3763,
                tmAlignAlignmentResult.getTemplateModelingScore1().getScore(),
                TestUtils.TOLERANT_ERROR_MARGIN);
        Assert.assertEquals("tmscore does not match",
                0.2059,
                tmAlignAlignmentResult.getTemplateModelingScore2().getScore(),
                TestUtils.TOLERANT_ERROR_MARGIN);
    }

    @Test(expected = AlignmentException.class)
    @Ignore
    public void shouldReportException() throws AlignmentException {
        instance.process(new String[] { "tmalign", "invalid", "arguments" });
    }
}
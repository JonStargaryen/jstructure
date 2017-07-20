package de.bioforscher.jstructure.align.impl;

import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.stream.Collectors;

/**
 * Test for the local blast wrapper.
 * Created by bittrich on 7/18/17.
 */
public class LocalBlastWrapperTest {
    private LocalBlastWrapper blastWrapper;
    private String sequence;

    @Before
    public void setup() {
        blastWrapper = new LocalBlastWrapper();
        sequence = TestUtils.getResourceAsStream("1acj.fsa")
                .collect(Collectors.joining(System.lineSeparator()));
    }

    @Test
    public void shouldExecutePsiBlast() {
        LocalBlastWrapper.PsiBlastResult result = blastWrapper.executePsiBlastUniref50(sequence);
        Assert.assertFalse(result.getAccessions().isEmpty());
        Assert.assertFalse(result.getInformation().isEmpty());
        Assert.assertEquals("number of conservation scores does not match",
                537,
                result.getInformation().size());
    }

    @Test
    public void shouldParseMatrixFile() {
        Assert.assertTrue("did not extract conservation scores",
                blastWrapper.parseMatrixFile(TestUtils.getResourceAsStream("1acj.mtx")).size() > 0);
    }

    @Test
    public void shouldParseOutputFile() {
        Assert.assertTrue("did not extract accessions",
                blastWrapper.parseResultFile(TestUtils.getResourceAsStream("1acj.out")).size() > 0);
    }
}
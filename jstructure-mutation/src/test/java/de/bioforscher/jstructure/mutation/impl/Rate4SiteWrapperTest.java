package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.stream.Collectors;

/**
 * Test the rate4site wrapper.
 * Created by bittrich on 7/17/17.
 */
public class Rate4SiteWrapperTest {
    private String alignmentString;
    private Rate4SiteWrapper rate4SiteWrapper;

    @Before
    public void setup() {
        alignmentString = TestUtils.getResourceAsStream("rate4site/seq.aln")
                .collect(Collectors.joining(System.lineSeparator()));
        rate4SiteWrapper = new Rate4SiteWrapper();
    }

    @Test
    public void shouldExecuteTask() throws IOException, InterruptedException {
        System.out.println(rate4SiteWrapper.executeCommand(alignmentString));
    }
}
package de.bioforscher.jstructure.si;

import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.util.stream.Collectors;

public class ConfoldServiceWorkerTest {
    private String sequence;
    private String ss;
    private String rr;

    @Before
    public void setup() {
        sequence = TestUtils.getResourceAsStream("confold/short.fasta")
                .collect(Collectors.joining(System.lineSeparator()));
        ss = TestUtils.getResourceAsStream("confold/short.ss")
                .collect(Collectors.joining(System.lineSeparator()));
        rr = TestUtils.getResourceAsStream("confold/short.rr")
                .collect(Collectors.joining(System.lineSeparator()));
    }

    @Test
    public void shouldRunReconstruction() {
        new ConfoldServiceWorker("/home/bittrich/programs/confold_v1.0/confold.pl",
                sequence,
                ss,
                rr)
                .call();
    }
}
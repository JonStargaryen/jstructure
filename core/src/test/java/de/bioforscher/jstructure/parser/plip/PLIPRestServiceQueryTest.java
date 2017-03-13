package de.bioforscher.jstructure.parser.plip;

import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;

/**
 * Test the PLIPRestServiceQuery.
 * Created by bittrich on 2/9/17.
 */
public class PLIPRestServiceQueryTest {
    @Test
    public void shouldReadCredentialsFile() {
        Assert.assertTrue(PLIPRestServiceQuery.secret != null);
    }

    @Test
    public void shouldRetrievePLIPDataSpecific() throws IOException {
        System.out.println(PLIPRestServiceQuery.getPlipResults(new URL(PLIPRestServiceQuery.BASE_URL + "1brr/A/100")));
    }

    @Test
    public void shouldRetrievePLIPDataForChain() throws IOException {
        System.out.println(PLIPRestServiceQuery.getPlipResults(new URL(PLIPRestServiceQuery.BASE_URL + "1brr/A")));
    }

    @Test
    public void shouldTestCoverageForDataSet() throws IOException {
        System.err.println("skipping PLIP-rest-service coverage test for data set");
//        Files.lines(Paths.get("/home/bittrich/git/phd_sb_repo/data/dataset/nrpdbtm/pdbtm_alpha_nr.list.txt"))
//                .map(line -> ProteinParser.source(line.substring(0, 4)).parse())
//                .forEach(PLIPRestServiceQuery::process);
    }
}
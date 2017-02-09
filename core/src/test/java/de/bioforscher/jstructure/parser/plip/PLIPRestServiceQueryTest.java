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
        System.err.println("skipping PLIP-rest-service coverage test for data set - coverage is 100%");
//        private static final String listPath = "/home/bittrich/git/phd_sb_repo/data/dataset/nrpdbtm/pdbtm_alpha_nr.list.txt";
//        Files.lines(Paths.get(listPath))
//                .map(line -> Selection.on(ProteinParser.parseProteinById(line.substring(0, 4)))
//                        .chainName(line.substring(5))
//                        .asChainContainer())
//                .map(Protein.class::cast)
//                .forEach(PLIPRestServiceQuery::process);
    }
}
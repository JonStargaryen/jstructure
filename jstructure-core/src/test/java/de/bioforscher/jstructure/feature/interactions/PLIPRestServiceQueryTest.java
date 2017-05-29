package de.bioforscher.jstructure.feature.interactions;

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
        System.out.println(PLIPRestServiceQuery.getDocument(new URL(PLIPRestServiceQuery.BASE_URL + "1brr/A/100")));
    }

    @Test
    public void shouldRetrievePLIPDataForChain() throws IOException {
        System.out.println(PLIPRestServiceQuery.getDocument(new URL(PLIPRestServiceQuery.BASE_URL + "1brr/A")));
    }
}
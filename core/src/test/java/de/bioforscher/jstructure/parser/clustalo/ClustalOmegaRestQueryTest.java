package de.bioforscher.jstructure.parser.clustalo;

import org.junit.Before;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Test the clustal omega REST interface.
 * Created by S on 17.03.2017.
 */
public class ClustalOmegaRestQueryTest {
    private List<String> sequences;

    @Before
    public void setup() {
        sequences = Stream.of("1M0L", "1AT9", "1BM1", "1BRD")
                .map(id -> {
                    try {
                        URL fetchUrl = new URL("http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=" + id);
                        try (BufferedReader reader = new BufferedReader(new InputStreamReader(fetchUrl.openStream(), StandardCharsets.UTF_8))) {
                            return reader.lines().collect(Collectors.joining(System.lineSeparator()));
                        }
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                })
                .collect(Collectors.toList());
    }

    @Test
    public void shouldPostClustalOmegaQuery() throws ExecutionException {
        //TODO test
        System.out.println(new ClustalOmegaRestQuery().process(sequences));
    }
}
package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.Base64;
import java.util.stream.Collectors;

/**
 * A wrapper for the PLIP-rest-service which provides precomputed PLIP results.
 * Created by bittrich on 2/9/17.
 */
public class PLIPRestServiceQuery {
    private static final Logger logger = LoggerFactory.getLogger(PLIPRestServiceQuery.class);
    static final String BASE_URL = "https://biosciences.hs-mittweida.de/plip/interaction/";
//    private static final String REST_USER_PASSWORD_PATH = System.getProperty("user.home") + "/git/phd_sb_repo/data/.plip-rest-auth";
    static String secret;

    static {
        try {
//            String line = Files.readAllLines(Paths.get(REST_USER_PASSWORD_PATH)).get(0);
            String line = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream("plip_credentials.txt"))).readLine();
            secret = new String(Base64.getMimeEncoder().encode(line.getBytes()));
        } catch (IOException e) {
            throw new IllegalStateException("no credentials provided to access 'biosciences.hs-mittweida.de/plip/'");
        }
    }

    static Document getDocument(String pdbId, String chainId, int residueNumber) {
        try {
            return getDocument(new URL(BASE_URL + pdbId + "/" + chainId + "/" + residueNumber));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Document getDocument(ChainIdentifier chainId) {
        try {
            return getDocument(new URL(BASE_URL + "plain/" + chainId.getProteinIdentifier().getPdbId() + "/" + chainId.getChainId()));
        } catch (IOException e) {
            throw new UncheckedIOException("failed to fetch PLIP files from REST service", e);
        }
    }

    public static Document getDocument(String pdbId, String chainId) {
        try {
            return getDocument(new URL(BASE_URL + "plain/" + pdbId + "/" + chainId));
        } catch (IOException e) {
            throw new UncheckedIOException("failed to fetch PLIP files from REST service", e);
        }
    }

    /**
     * Query the PLIP-rest-service for a given URL (which contains the concrete query).
     * @param url the query - currently supported are /pdbid/CHAINID/ and /pdbid/CHAINID/resnum
     * @return the server's answer as String
     * @throws IOException thrown when the service cannot be reached, authentication fails or the query is not supported
     */
    static Document getDocument(URL url) throws IOException {
        logger.info("querying PLIP-rest-service for {}", url.toString());
        URLConnection connection = url.openConnection();
        connection.setRequestProperty("Authorization", "Basic " + secret);
        connection.connect();

        try (InputStream inputStream = connection.getInputStream()) {
            try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream))) {
                return Jsoup.parse(bufferedReader.lines()
                        .collect(Collectors.joining(System.lineSeparator())));
            }
        }
    }
}

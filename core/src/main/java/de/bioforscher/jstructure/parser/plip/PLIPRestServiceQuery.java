package de.bioforscher.jstructure.parser.plip;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

/**
 * A wrapper for the PLIP-rest-service which provides precomputed PLIP results.
 * Created by bittrich on 2/9/17.
 */
class PLIPRestServiceQuery {
    private static final Logger logger = LoggerFactory.getLogger(PLIPRestServiceQuery.class);
    static final String BASE_URL = "http://141.55.231.200:8731/plip/";
    private static final String REST_USER_PASSWORD_PATH = "/home/bittrich/git/phd_sb_repo/data/.plip-rest-auth";
    static String secret;
    private static final PLIPRestServiceQuery INSTANCE = new PLIPRestServiceQuery();

    private PLIPRestServiceQuery() {
        try {
            String line = Files.readAllLines(Paths.get(REST_USER_PASSWORD_PATH)).get(0);
            secret = new sun.misc.BASE64Encoder().encode(line.getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static String getPlipResults(String pdbId, String chainId, int residueNumber) {
        try {
            return getPlipResults(new URL(BASE_URL + pdbId + "/" + chainId + "/" + residueNumber));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static String getPlipResults(String pdbId, String chainId) {
        try {
            return getPlipResults(new URL(BASE_URL + pdbId + "/" + chainId));
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
    static String getPlipResults(URL url) throws IOException {
        logger.info("querying PLIP-rest-service for {}", url.toString());
        URLConnection connection = url.openConnection();
        connection.setRequestProperty("Authorization", "Basic " + secret);
        connection.connect();

        try (InputStream inputStream = connection.getInputStream()) {
            try (BufferedReader buffer = new BufferedReader(new InputStreamReader(inputStream))) {
                return buffer.lines().collect(Collectors.joining(System.lineSeparator()));
            }
        }
    }
}
package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Base64;
import java.util.stream.Collectors;

/**
 * A wrapper for the PLIP-rest-service which provides precomputed PLIP results.
 * Created by bittrich on 2/9/17.
 */
public class PLIPRestServiceQuery {
    private static final Logger logger = LoggerFactory.getLogger(PLIPRestServiceQuery.class);
    static final String BASE_URL = "https://biosciences.hs-mittweida.de/plip/interaction/plain/";
    static String secret;

    static {
        try {
            String line = new BufferedReader(new InputStreamReader(Thread.currentThread()
                    .getContextClassLoader()
                    .getResourceAsStream("plip_credentials.txt"))).readLine();
            secret = new String(Base64.getMimeEncoder().encode(line.getBytes()));
            logger.info("PLIP Service is running against: {} - authentication provided", BASE_URL);
        } catch (IOException | NullPointerException e) {
            throw new IllegalStateException("no credentials provided to access 'biosciences.hs-mittweida.de/plip/'");
        }
    }

    public static Document getIntraMolecularDocument(Chain chain) {
        ChainIdentifier chainIdentifier = chain.getChainIdentifier();
        return getIntraMolecularDocument(chainIdentifier.getProteinIdentifier().getPdbId(),
                chainIdentifier.getChainId());
    }

    /**
     * Compute interactions for a chain not previously processed. Designed for CASP data.
     * @param chain the chain to process
     * @return the document containing all results
     */
    public static Document calculateIntraChainDocument(Chain chain) {
        try {
            String chainId = chain.getChainIdentifier().getChainId();
            // write PDB structure of data point to temporary file
            Path structureFilePath = Files.createTempFile("plip_", "_" + chainId + ".pdb");
            System.out.println(chain.getPdbRepresentation());
            Files.write(structureFilePath, chain.getPdbRepresentation().getBytes());

            // submit PLIP POST query
            String url = "https://biosciences.hs-mittweida.de/plip/interaction/calculate/intrachain/" + chainId;
            logger.info("processing chain '{}' by remote PLIP at:{}'{}'",
                    chainId,
                    System.lineSeparator(),
                    url);
            PLIPPostRequest plipPostRequest = new PLIPPostRequest(url,
                    secret,
                    structureFilePath);
            return Jsoup.parse(plipPostRequest.getResultXml());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static Document getLigandDocument(Chain chain) {
        try {
            // write PDB structure of data point to temporary file
            Path structureFilePath = Files.createTempFile("plip_", "_" + chain.getChainIdentifier().getFullName() + ".pdb");
            Files.write(structureFilePath, chain.getPdbRepresentation().getBytes());

            // submit PLIP POST query
            PLIPPostRequest plipPostRequest = new PLIPPostRequest("https://biosciences.hs-mittweida.de/plip/interaction/calculate/protein", secret, structureFilePath);
            return Jsoup.parse(plipPostRequest.getResultXml());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static Document getIntraMolecularDocument(String pdbId, String chainId) {
        try {
            return getIntraMolecularDocument(new URL(BASE_URL + pdbId + "/" + chainId));
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
    static Document getIntraMolecularDocument(URL url) throws IOException {
        logger.debug("querying PLIP-rest-service for {}", url.toString());
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

    static class PLIPPostRequest {
//        private static final Logger logger = LoggerFactory.getLogger(PLIPPostRequest.class);
        private final String credentials;
        private final Path pdbFilePath;
        private String plipUrl;
        private String resultXml;

        PLIPPostRequest(String plipUrl, String credentials, Path pdbFilePath) throws IOException {
            this.plipUrl = plipUrl;
            this.credentials = credentials;
            this.pdbFilePath = pdbFilePath;
//            logger.info("creating PLIP POST for PDB file {}", pdbFilePath);
            this.doPlipPost();
        }

        String getResultXml() {
            return this.resultXml;
        }

        private void doPlipPost() throws IOException {
            String boundary = Long.toHexString(System.currentTimeMillis());
            String crlf = "\r\n";
            URLConnection connection = (new URL(this.plipUrl)).openConnection();
            connection.setDoOutput(true);
            connection.setRequestProperty("Content-Type", "multipart/form-data; boundary=" + boundary);
            connection.setRequestProperty("Authorization", "Basic " + credentials);
            OutputStream output = connection.getOutputStream();
            PrintWriter writer = new PrintWriter(new OutputStreamWriter(output, "UTF-8"), true);
            writer.append("--").append(boundary).append(crlf);
            writer.append("Content-Disposition: form-data; name=\"file\"; filename=\"").append(this.pdbFilePath.getFileName().toString()).append("\"").append(crlf);
            writer.append("Content-Type: text/plain; charset=UTF-8").append(crlf);
            writer.append(crlf).flush();
            Files.copy(this.pdbFilePath, output);
            output.flush();
            writer.append(crlf).flush();
            writer.append("--").append(boundary).append("--").append(crlf).flush();
//            int responseCode = ((HttpURLConnection)connection).getResponseCode();
//            logger.info("PLIP REST service response code is " + responseCode);
            BufferedReader in = new BufferedReader(new InputStreamReader(connection.getInputStream()));
            StringBuilder response = new StringBuilder();

            String outputContent;
            while((outputContent = in.readLine()) != null) {
                response.append(outputContent);
            }

            in.close();
            this.resultXml = response.toString();
        }
    }
}

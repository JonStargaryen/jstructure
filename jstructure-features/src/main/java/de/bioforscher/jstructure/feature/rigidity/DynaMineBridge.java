package de.bioforscher.jstructure.feature.rigidity;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * see: http://dynamine.ibsquare.be/download/
 */
public class DynaMineBridge extends FeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(DynaMineBridge.class);
    private final String apiKey;
    private static final String BASE_URL = "http://dynamine.ibsquare.be/batch_request";

    public DynaMineBridge() {
        try {
            this.apiKey = new BufferedReader(new InputStreamReader(Thread.currentThread()
                    .getContextClassLoader()
                    .getResourceAsStream("dynamine_credentials.txt"))).readLine();
            logger.info("DynaMine Service is running against: {} - authentication provided",
                    BASE_URL);
        } catch (IOException e) {
            logger.warn("loading DynaMine credentials file failed - provide your own api-key as resource named 'dynamine_credentials.txt'");
            throw new UncheckedIOException(e);
        }
    }

    @Override
    protected void processInternally(Structure structure) {
        structure.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        String sequence = chain.getAminoAcidSequence();
        logger.info("processing chain {} with sequence {}",
                chain.getChainIdentifier(),
                sequence);

        try {
            DynaMineRequest request = new DynaMineRequest(sequence);
            String payload = request.asJson();
            logger.info("payload is:{}{}",
                    System.lineSeparator(),
                    payload);
            byte[] out = payload.getBytes(StandardCharsets.UTF_8);
            int length = out.length;

            URL url = new URL(BASE_URL);
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("POST");
            connection.setDoOutput(true);

            connection.setFixedLengthStreamingMode(length);
            connection.setRequestProperty("Content-Type", "application/json; charset=UTF-8");
            connection.connect();

            try(OutputStream os = connection.getOutputStream()) {
                os.write(out);
            }

            int responseCode = connection.getResponseCode();
            if(responseCode != 200) {
                logger.warn("response code was {}{}{}",
                        responseCode,
                        System.lineSeparator(),
                        new BufferedReader(new InputStreamReader(connection.getErrorStream()))
                        .lines()
                        .collect(Collectors.joining(System.lineSeparator())));
            }

            //TODO server returns 500 - can this be fixed?
            new BufferedReader(new InputStreamReader(connection.getInputStream()))
                    .lines()
                    .forEach(System.out::println);
        } catch (Exception e) {
            logger.warn("submitting job to DynaMine failed",
                    e);
        }
    }

    public void process(Chain chain, String dynamineDocument) {
        List<Double> predictions =  Pattern.compile("[\\r\\n]+").splitAsStream(dynamineDocument)
                .filter(line -> line.contains("\t"))
                .map(line -> line.split("\t")[1])
                .map(Double::valueOf)
                .collect(Collectors.toList());
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

        double min = predictions.stream()
                .mapToDouble(Double::valueOf)
                .min()
                .orElse(0);
        double max = predictions.stream()
                .mapToDouble(Double::valueOf)
                .max()
                .orElse(1);

        if(predictions.size() != aminoAcids.size()) {
            throw new IllegalArgumentException("parsed lines do not match expectation - " + predictions.size() + " for " +
                    aminoAcids.size() + " amino acids");
        }

        for (int i = 0; i < predictions.size(); i++) {
            AminoAcid aminoAcid = aminoAcids.get(i);
            double value = minMaxNormalize(predictions.get(i), min, max);
            aminoAcid.getFeatureContainer().addFeature(new BackboneRigidity(this, value));
        }
    }

    private static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }

    class DynaMineRequest {
        private final String protocol = "1.0";
        private final String json_api_key;
        private final String sequence;
        // false will produce additional plots etc
        private final boolean predictions_only = true;

        DynaMineRequest(String sequence) {
            this.json_api_key = DynaMineBridge.this.apiKey;
            this.sequence = sequence;
        }

        String asJson() {
            return "{" + System.lineSeparator() +
                    "\"protocol\" : \"" + protocol + "\"," + System.lineSeparator() +
                    "\"json_api_key\" : \"" + json_api_key + "\"," + System.lineSeparator() +
                    "\"sequences\" : { \"P04637\" : \"" + sequence + "\" }," + System.lineSeparator() +
                    "\"predictions_only\" : " + predictions_only + System.lineSeparator() +
                    "}";
        }
    }
}

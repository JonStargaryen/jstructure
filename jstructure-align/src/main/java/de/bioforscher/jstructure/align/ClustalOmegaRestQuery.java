package de.bioforscher.jstructure.align;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Queries the REST interface of ebi to run clustal omega.
 * TODO wrap in services
 * Created by bittrich on 7/11/17.
 */
public class ClustalOmegaRestQuery {
    private static final Logger logger = LoggerFactory.getLogger(ClustalOmegaRestQuery.class);
    private static final String BASE_URL = "http://www.ebi.ac.uk/Tools/services/rest/clustalo/";
    private static final String RUN_URL = BASE_URL + "run/";
    private static final String STATUS_URL = BASE_URL + "status/%s";
    private static final String RESULT_URL = BASE_URL + "result/%s/aln-fasta";
    private static final Pattern SEQUENCE_PATTERN = Pattern.compile(">");
    private static final Pattern LINE_PATTERN = Pattern.compile("\\s+");

    public ClustalOmegaRestQuery() {
        logger.info("ClustalOmega Service is running against {}", RUN_URL);
    }

    /**
     * Computes a multiple-sequence alignment for a given collection of sequences.
     * @param sequences the input sequences - expected to be protein sequence in FASTA format
     * @return a map with the id as key and the aligned sequence as value
     * @throws ExecutionException upon failed execution
     */
    public Map<String, String> process(List<String> sequences) throws ExecutionException {
        sequences.forEach(sequence -> {
            if(!sequence.contains(">") || LINE_PATTERN.matcher(sequence).find()) {
                throw new IllegalArgumentException("sequence must be in FASTA format - found: '" + sequence + "'");
            }
        });

        logger.info("creating multi-sequence alignment by clustal omega for {} protein sequences", sequences.size());
        try {
            Document answer = Jsoup.connect(RUN_URL)
                    .data("sequence", composeSequenceString(sequences))
                    .ignoreHttpErrors(true)
                    //TODO remove
                    .data("email", "bittrich@hs-mittweida.de")
                    .data("outfmt", "fa")
                    .data("stype", "protein")
                    .post();

            String jobId = answer.text();
            logger.info("job id is {}", jobId);
            waitForResults(jobId);

            String rawAlignment = Jsoup.connect(String.format(RESULT_URL, jobId))
                    .get()
                    .text();

            Map<String, String> alignment = new HashMap<>();
            SEQUENCE_PATTERN.splitAsStream(rawAlignment)
                    // skip first (empty) pair
                    .skip(1)
                    .forEach(line -> {
                        // split alignment at newlines
                        String[] split = LINE_PATTERN.split(line);
                        // substring to drop FASTA-character
                        alignment.put(split[0].substring(1),
                                // skip id and join all other lines
                                Stream.of(split).skip(1).collect(Collectors.joining()));
                    });

            return alignment;
        } catch (InterruptedException | IOException | RuntimeException e) {
            throw new ExecutionException(e);
        }
    }

    private void waitForResults(String jobId) throws IOException, InterruptedException {
        String status = Jsoup.connect(String.format(STATUS_URL, jobId)).get().text();

        // error case
        if(!status.equals("RUNNING") && !status.equals("FINISHED")) {
            //TODO consistent error handling
            throw new RuntimeException("job did not run or finish");
        }

        if(!status.equals("FINISHED")) {
            TimeUnit.SECONDS.sleep(1);
            waitForResults(jobId);
        }
    }

    private String composeSequenceString(List<String> sequences) {
        return sequences.stream()
                .collect(Collectors.joining(System.lineSeparator()));
    }
}

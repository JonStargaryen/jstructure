package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Queries the REST interface of ebi to run clustal omega.
 * Created by bittrich on 7/11/17.
 */
public class ClustalOmegaRestQuery implements MultipleSequenceAligner {
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

    @Override
    public MultipleSequenceAlignmentResult align(List<String> fastaSequences) throws AlignmentException {
        fastaSequences.forEach(sequence -> {
            if(!sequence.startsWith(">")) {
                throw new IllegalArgumentException("sequence must be in FASTA format - found:" + System.lineSeparator()
                        + sequence);
            }
        });

        logger.info("creating multi-sequence alignment by ClustalOmega for {} protein sequences", fastaSequences.size());
        try {
            Document answer = Jsoup.connect(RUN_URL)
                    .data("sequence", composeSequenceString(fastaSequences))
                    .ignoreHttpErrors(true)
                    //TODO remove
                    .data("email", "bittrich@hs-mittweida.de")
                    .data("outfmt", "fa")
                    .data("stype", "protein")
                    .post();

            String jobId = answer.text();
            logger.debug("job id is {}", jobId);
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
                        alignment.put(split[0],
                                // skip id and join all other lines
                                Stream.of(split).skip(1).collect(Collectors.joining()));
                    });

            return new MultipleSequenceAlignmentResultImpl(alignment);
        } catch (InterruptedException | IOException | RuntimeException e) {
            throw new AlignmentException(e);
        }
    }

    private void waitForResults(String jobId) throws IOException, InterruptedException {
        String status = Jsoup.connect(String.format(STATUS_URL, jobId)).get().text();

        // error case
        if(!status.equals("RUNNING") && !status.equals("FINISHED")) {
            throw new AlignmentException("job did not run or finish");
        }

        if(!status.equals("FINISHED")) {
            TimeUnit.SECONDS.sleep(5);
            waitForResults(jobId);
        }
    }

    private String composeSequenceString(List<String> sequences) {
        return sequences.stream()
                .collect(Collectors.joining(System.lineSeparator()));
    }
}

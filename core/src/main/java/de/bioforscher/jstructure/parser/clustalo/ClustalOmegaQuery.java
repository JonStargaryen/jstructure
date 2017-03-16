package de.bioforscher.jstructure.parser.clustalo;

import com.gargoylesoftware.htmlunit.WebClient;
import com.gargoylesoftware.htmlunit.html.HtmlOption;
import com.gargoylesoftware.htmlunit.html.HtmlPage;
import com.gargoylesoftware.htmlunit.html.HtmlSelect;
import com.gargoylesoftware.htmlunit.html.HtmlTextArea;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * Posts a MSA-query to clustal omega.
 * Created by bittrich on 3/16/17.
 */
public class ClustalOmegaQuery {
    private static final Logger logger = LoggerFactory.getLogger(ClustalOmegaQuery.class);
    public static final String CLUSTAL_OMEGA_ALIGNMENT = "CLUSTAL_OMEGA_ALIGNMENT";
    private static final String BASE_URL = "http://www.ebi.ac.uk/Tools/msa/clustalo/";
    private static final String RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/clustalo/result/%s/aln-fasta";
    private static final int WAIT_SECONDS = 5;

    public String process(List<String> sequences) {
        logger.info("creating multi-sequence alignment by clustal omega for {} proteins", sequences.size());

        try(final WebClient webClient = new WebClient()) {
            try {
                final HtmlPage form = webClient.getPage(BASE_URL);
                ((HtmlTextArea) form.getElementById("sequence")).setText(composeSequenceString(sequences));

                HtmlSelect select = (HtmlSelect) form.getElementById("outfmt");
                HtmlOption option = select.getOptionByValue("fa");
                select.setSelectedAttribute(option, true);

                final HtmlPage queue = form.getElementById("jd_submitButtonPanel").getFirstElementChild().click();
                String jobId = queue.asText().split("Job ID: ")[1].split("[\\r\\n]+")[0];
                logger.debug("jobid is {}", jobId);

                URL resultUrl = new URL(String.format(RESULT_URL, jobId));
                waitForResults(resultUrl);

                URLConnection conn = resultUrl.openConnection();
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream(), StandardCharsets.UTF_8))) {
                    String resultText = reader.lines().collect(Collectors.joining(System.lineSeparator()));
                    logger.debug("alignment result:{}{}", System.lineSeparator(), resultText);
                    return resultText;
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }

    private void waitForResults(URL url) {
        try {
            url.getContent();
        } catch (IOException e) {
            logger.trace("results not available yet");
            try {
                TimeUnit.SECONDS.sleep(WAIT_SECONDS);
            } catch (InterruptedException e1) {
                logger.error("interrupted while waiting for clustal omega results");
            }
            waitForResults(url);
        }
    }

    String composeSequenceString(List<String> sequences) {
        return sequences.stream()
                .collect(Collectors.joining(System.lineSeparator()));
    }

//    String composeSequenceString(List<Protein> proteins) {
//        return proteins.stream()
//                .flatMap(Protein::chains)
//                // keep only chains with amino acids
//                .filter(chain -> chain.aminoAcids().count() > 0)
//                .map(chain -> ">" + chain.getParentProtein().getName() + "." + chain.getChainId() + System.lineSeparator() + chain.getAminoAcidSequence())
//                .collect(Collectors.joining(System.lineSeparator()));
//    }
}

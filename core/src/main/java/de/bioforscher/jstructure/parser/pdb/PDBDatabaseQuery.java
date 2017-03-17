package de.bioforscher.jstructure.parser.pdb;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Access to pdb functionality provided via REST.
 * Created by bittrich on 3/17/17.
 */
public class PDBDatabaseQuery {
    private static final Logger logger = LoggerFactory.getLogger(PDBDatabaseQuery.class);
    public static final String CLUSTER_FETCH_URL = "http://www.rcsb.org/pdb/rest/sequenceCluster?cluster=%d&structureId=%s.%s";
    public static final double DEFAULT_SEQUENCE_SIMILARITY_THRESHOLD = 95;

    /**
     * @see PDBDatabaseQuery#fetchSequenceCluster(String, String, double)
     */
    public static List<String> fetchSequenceCluster(String pdbId, String chainId) {
        return fetchSequenceCluster(pdbId, chainId, DEFAULT_SEQUENCE_SIMILARITY_THRESHOLD);
    }

    /**
     * Queries the PDB for all chains which share a given sequence similarity to a reference chain.
     * @param pdbId the reference pdb id
     * @param chainId the reference chain id
     * @param sequenceSimilarityThreshold the desired sequence similary threshold
     * @return all chain ids (in the format pdbId.chainId) which are similar to the reference
     */
    public static List<String> fetchSequenceCluster(String pdbId, String chainId, double sequenceSimilarityThreshold) {
        try {
            Document document = Jsoup.connect(String.format(CLUSTER_FETCH_URL, (int) sequenceSimilarityThreshold, pdbId, chainId)).get();
            return document.getElementsByTag("pdbChain").stream()
                    .map(pdbChain -> pdbChain.attr("name"))
                    .collect(Collectors.toList());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}

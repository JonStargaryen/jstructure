package de.bioforscher.explorer.membrane.model.homologous;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.clustalo.ClustalOmegaQuery;
import de.bioforscher.jstructure.parser.opm.OPMDatabaseQuery;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Information on homologous chains for a given reference protein.
 * Created by bittrich on 3/15/17.
 */
public class HomologousChainContainer {
    private static final Logger logger = LoggerFactory.getLogger(HomologousChainContainer.class);
    private List<HomologousChain> alignedChains;
    private int[][] comparison;
    private static final char GAP_CHAR = '-';
    public static final int GAP_VALUE = 0;
    public static final int MISMATCHING_GAP_VALUE = 4;
    public static final int MATCHING_GAP_VALUE = 8;
    public static final int MISMATCH_VALUE = 12;
    public static final int MATCH_VALUE = 16;

    public HomologousChainContainer() {
    }

    public HomologousChainContainer(Protein protein) {
        List<String> ids = new ArrayList<>();
        String pdbId = protein.getName();
        ids.add(pdbId);
        protein.getFeatureAsList(String.class, OPMDatabaseQuery.HOMOLOGOUS_PROTEINS).forEach(ids::add);

        logger.info("sequence cluster for {} is {}", pdbId, ids);

        List<String> homologousSequences = ids.stream()
                .map(this::fetchFastaSequence)
                .collect(Collectors.toList());

        this.alignedChains = convert(new ClustalOmegaQuery().process(homologousSequences), pdbId);

        this.comparison = computeComparisonArray(alignedChains);
    }

    private int[][] computeComparisonArray(List<HomologousChain> alignedChains) {
        HomologousChain referenceChain = alignedChains.get(0);
        String referenceSequence = referenceChain.getSequence();
        int[][] comparisonArray = new int[alignedChains.size()][];

        comparisonArray[0] = new int[referenceSequence.length()];
        for(int position = 0; position < referenceSequence.length(); position++) {
            comparisonArray[0][position] = referenceSequence.charAt(position) == GAP_CHAR ? GAP_VALUE : MATCH_VALUE;
        }

        for(int sequenceIndex = 1; sequenceIndex < alignedChains.size(); sequenceIndex++) {
            comparisonArray[sequenceIndex] = computeComparisonArray(referenceSequence, alignedChains.get(sequenceIndex).getSequence());
        }

        return comparisonArray;
    }

    private int[] computeComparisonArray(String referenceSequence, String homologousSequence) {
        int[] comparisonArray = new int[referenceSequence.length()];
        for(int position = 0; position < referenceSequence.length(); position++) {
            char referenceChar = referenceSequence.charAt(position);
            char homologousChar = homologousSequence.charAt(position);

            if(referenceChar == GAP_CHAR && homologousChar == GAP_CHAR) {
                // both gaps
                comparisonArray[position] = MATCHING_GAP_VALUE;
            } else if(referenceChar == GAP_CHAR || homologousChar == GAP_CHAR) {
                // one sequence gap
                comparisonArray[position] = MISMATCHING_GAP_VALUE;
            } else if(referenceChar == homologousChar) {
                // both sequence identical amino acid
                comparisonArray[position] = MATCH_VALUE;
            } else {
                comparisonArray[position] = MISMATCH_VALUE;
            }
        }
        return comparisonArray;
    }

    private List<HomologousChain> convert(String alignmentString, String pdbId) {
        List<String> split = Pattern.compile(">")
                .splitAsStream(alignmentString)
                .skip(1)
                .collect(Collectors.toList());

        List<HomologousChain> sequences = split.stream()
                .map(s -> s.split("SEQUENCE"))
                // substring for pdbId:chainId
                .map(s -> new Pair<>(s[0].substring(0, 6), Pattern.compile("[\\r\\n]+").splitAsStream(s[1]).collect(Collectors.joining())))
                .map(pair -> new HomologousChain(pair.getLeft(), pair.getRight()))
                .peek(System.out::println)
                .collect(Collectors.toList());

        // move representative sequence to the beginning - can there by more than 1?
        List<HomologousChain> representatives = sequences.stream()
                .filter(chain -> chain.getId().startsWith(pdbId))
                .collect(Collectors.toList());
        sequences.removeAll(representatives);
        representatives.forEach(chain -> sequences.add(0, chain));

        return sequences;
    }

    private String fetchFastaSequence(String pdbId) {
        try {
            URL fetchUrl = new URL("http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=" + pdbId);
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(fetchUrl.openStream(), StandardCharsets.UTF_8))) {
                return reader.lines().collect(Collectors.joining(System.lineSeparator()));
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public List<HomologousChain> getChains() {
        return alignedChains;
    }

    public int[][] getComparison() {
        return comparison;
    }
}

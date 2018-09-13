package de.bioforscher.jstructure.efr.analysis;

import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A01_SingleLinkageClusterer {
    private static final Logger logger = LoggerFactory.getLogger(A01_SingleLinkageClusterer.class);
    public static final double SIMILARITY_THRESHOLD = 0.95;
    private List<Chain> chains;
    private List<List<Chain>> clusters;

    A01_SingleLinkageClusterer() throws IOException {
        this.chains = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(line -> {
                    String[] split = line.split(";");
                    String entryId = split[0];
                    String pdbId = split[1];
                    Structure structure = StructureParser.fromPdbId(pdbId).parse();
                    return structure.chains().findFirst().get();
                })
                .collect(Collectors.toList());
        this.clusters = new ArrayList<>();
        this.cluster();

        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        for (int i = 0; i < clusters.size(); i++) {
            List<Chain> clusterContent = clusters.get(i);
            int clusterIndex = (i + 1);
            String clusterString = clusterIndex + " : " + clusterContent.stream()
                    .map(Chain::getChainIdentifier)
                    .map(ChainIdentifier::toString)
                    .collect(Collectors.toList());
            logger.info(clusterString);
            stringJoiner.add(clusterString);
        }
        System.out.println(stringJoiner.toString());
    }

    public static void main(String[] args) throws IOException {
        new A01_SingleLinkageClusterer();
    }

    private void cluster() {
        chains.forEach(sequenceToCheck -> {
            // stores indices of clusters which contain highly similar sequences
            List<Integer> indicesOfSimilarClusters = new ArrayList<>();

            for (int sequenceClusterIndex = 0; sequenceClusterIndex < clusters.size(); sequenceClusterIndex++) {
                boolean clusterContainsSimilarSequence = false;
                List<Chain> sequenceCluster = clusters.get(sequenceClusterIndex);
                for (Chain clusterEntry : sequenceCluster) {
                    // sequences are highly similar
                    if (computeNeedlemanWunschSimilarity(sequenceToCheck, clusterEntry) > SIMILARITY_THRESHOLD) {
                        clusterContainsSimilarSequence = true;
                        break;
                    }
                }

                if (clusterContainsSimilarSequence) {
                    indicesOfSimilarClusters.add(sequenceClusterIndex);
                }
            }

            if (indicesOfSimilarClusters.isEmpty()) {
                logger.info("sequence is unique like you, creating a new cluster {}", (clusters.size() + 1));
                // if no similar sequence clusters were found, we create a new cluster
                clusters.add(Stream.of(sequenceToCheck).collect(Collectors.toList()));
            } else {
                // the hard case: add to single cluster or merge clusters
                if (indicesOfSimilarClusters.size() == 1) {
                    logger.info("added sequence to cluster: " + (indicesOfSimilarClusters.get(0) + 1));
                    // the easy, yet hard case: add sequence to existing cluster
                    clusters.get(indicesOfSimilarClusters.get(0)).add(sequenceToCheck);
                } else {
                    // we need to join/merge existing clusters, because the processed sequence 'connects' them
                    logger.info("merging clusters with indices: " + indicesOfSimilarClusters);

                    // add first (a.k.a. new/currently processed sequence)
                    List<Chain> mergedCluster = Stream.of(sequenceToCheck).collect(Collectors.toList());
                    // zeck the similar clusters and add them to the joint cluster
                    for (int indexOfSimilarCluster : indicesOfSimilarClusters) {
                        mergedCluster.addAll(clusters.get(indexOfSimilarCluster));
                    }
                    // maybe we do not have to sort, as indices should be in ascending order, but just to be sure
                    Collections.sort(indicesOfSimilarClusters);
                    Collections.reverse(indicesOfSimilarClusters);
                    for (int i : indicesOfSimilarClusters) {
                        clusters.remove(i);
                    }
                }
            }
        });
    }

    /**
     * Compute the Needleman-Wunsch alignment between 2 sequence and report the sequence identity.
     *
     * @param entry1 the reference sequence
     * @param entry2 the query sequence
     * @return the fraction of identically aligned positions
     */
    private double computeNeedlemanWunschSimilarity(Chain entry1, Chain entry2) {
        try {
            ProteinSequence sequence1 = new ProteinSequence(entry1.getAminoAcidSequence());
            ProteinSequence sequence2 = new ProteinSequence(entry2.getAminoAcidSequence());

            SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(sequence1,
                    sequence2,
                    Alignments.PairwiseSequenceAlignerType.GLOBAL,
                    new SimpleGapPenalty(),
                    SubstitutionMatrixHelper.getBlosum62());

            System.out.println(pair.getPercentageOfIdentity());
            return pair.getPercentageOfIdentity();
        } catch (CompoundNotFoundException e) {
            throw new IllegalArgumentException(e);
        }
    }
}

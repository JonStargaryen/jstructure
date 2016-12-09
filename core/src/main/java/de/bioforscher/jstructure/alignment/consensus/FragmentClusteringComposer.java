package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * An clustering approach somewhat diametric to {@link ConsensusTreeComposer}.
 * This time, fragments are assigned to clusters when they are highly similar to it and move to novel ones, when they do
 * not share significant similarity to an existing cluster.
 * The consensus fragment is composed by merging the existing consensus fragment with that newly assigned to this
 * cluster.
 * Created by S on 05.12.2016.
 */
public class FragmentClusteringComposer extends AbstractConsensusComposer {
    private static final double DEFAULT_RMSD_THRESHOLD = 1.0;
    private final double rmsdThreshold;
    private List<StructureCluster> clusters;

    public FragmentClusteringComposer() {
        this(DEFAULT_RMSD_THRESHOLD);
    }

    public FragmentClusteringComposer(double rmsdThreshold) {
        this.clusters = new ArrayList<>();
        this.rmsdThreshold = rmsdThreshold;
    }

    public void composeClusterRepresentation(List<? extends AtomContainer> containers) {
        for(AtomContainer container : containers) {
            try {
                List<StructureCluster> similarClusters = clusters.parallelStream()
                        .filter(cluster -> cluster.getDistanceToConsensus(container) < rmsdThreshold)
                        .collect(Collectors.toList());

                if (similarClusters.size() == 0) {
                    // not similar to any cluster
                    this.clusters.add(new StructureCluster(container));
                } else if (similarClusters.size() == 1) {
                    // similar to a single cluster
                    similarClusters.get(0).add(container);
                } else {
                    // remove merged clusters from collection
                    clusters.removeAll(similarClusters);

                    // similar to multiple clusters - move everything to first cluster
                    StructureCluster mergedCluster = similarClusters.remove(0);
                    // assign each other
                    similarClusters.stream()
                            .map(StructureCluster::getEntries)
                            .flatMap(Collection::stream)
                            .forEach(mergedCluster::add);

                    clusters.add(mergedCluster);
                }
            } catch (IllegalArgumentException e) {
                e.printStackTrace();
            }
        }

        System.out.println("grouped " + clusters.stream()
                .map(StructureCluster::getEntries)
                .mapToInt(Collection::size)
                .sum() + " in " + clusters.size() + " clusters");
    }

    public List<StructureCluster> getClusters() {
        return clusters;
    }

    public class StructureCluster {
        private List<AtomContainer> entries;
        private AtomContainer consensus;
        private SVDSuperimposer svdSuperimposer;

        StructureCluster(AtomContainer container) {
            this.entries = new ArrayList<>();
            this.entries.add(Objects.requireNonNull(container));
            this.consensus = container;
            this.svdSuperimposer = new SVDSuperimposer();
        }

        AtomContainer add(AtomContainer container) {
            entries.add(container);
            consensus = computeConsensus();
            return consensus;
        }

        public AtomContainer getConsensusRepresentation() {
            return consensus;
        }

        public List<AtomContainer> getEntries() {
            return entries;
        }

        private AtomContainer computeConsensus() {
            AlignmentResult alignment = svdSuperimposer.align(consensus, entries.get(entries.size() - 1));
            // merge new entry to existing consensus
            //TODO this may be fast, but potentially merging everything is more robust
            return mergeContainerPair(alignment.getOriginalReference(), alignment.getOriginalQuery());
        }

        double getMinimalDistance(AtomContainer container) {
            return entries.stream()
                    .map(entry -> svdSuperimposer.align(entry, container))
                    .mapToDouble(AlignmentResult::getRmsd)
                    .min()
                    .orElseThrow(() -> new IllegalArgumentException("rmsd to clusters could not be computed"));
        }

        double getDistanceToConsensus(AtomContainer container) {
            return svdSuperimposer.align(consensus, container).getRmsd();
        }
    }
}

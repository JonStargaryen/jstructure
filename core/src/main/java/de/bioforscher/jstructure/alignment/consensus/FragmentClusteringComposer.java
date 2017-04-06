package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An clustering approach somewhat diametric to {@link ConsensusTreeComposer}.
 * This time, fragments are assigned to clusters when they are highly similar to it and move to novel ones, when they do
 * not share significant similarity to an existing cluster.
 * The consensus fragment is composed by merging the existing consensus fragment with that newly assigned to this
 * cluster.
 * Created by S on 05.12.2016.
 */
@Deprecated
public class FragmentClusteringComposer extends AbstractConsensusComposer {
    private static final Logger logger = LoggerFactory.getLogger(FragmentClusteringComposer.class);
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
        for(int containerIndex = 0; containerIndex < containers.size(); containerIndex++) {
            AtomContainer container = containers.get(containerIndex);
            logger.info("handling container {} of {} with id: {}", containerIndex, containers.size(), container.getIdentifier());
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
                            .map(StructureCluster::getOriginalEntries)
                            .flatMap(Collection::stream)
                            .forEach(mergedCluster::add);

                    // 12/14/16 - fix: forgot to add current container to merged containers
                    mergedCluster.add(container);

                    clusters.add(mergedCluster);
                }
            } catch (IllegalArgumentException e) {
                e.printStackTrace();
            }
        }

        logger.info("grouped {} in {} clusters", clusters.stream()
                .map(StructureCluster::getOriginalEntries)
                .mapToInt(Collection::size)
                .sum(), clusters.size());
    }

    public List<StructureCluster> getClusters() {
        return clusters;
    }
}

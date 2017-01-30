package de.bioforscher.jstructure.alignment.multiple;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.alignment.consensus.AbstractConsensusComposer;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Performs a multiple structure alignment.
 * Created by bittrich on 1/25/17.
 */
public class MultipleStructureAlignment {
    private static final Logger logger = LoggerFactory.getLogger(MultipleStructureAlignment.class);
    private final SVDSuperimposer svdSuperimposer;
    /**
     * The cutoff when a binding site is considered to not be able to extend any further.
     */
    private static final double RMSD_CUTOFF = 2.0;
    private static final int NUMBER_OF_GROUPS_AS_CANDIDATES_FOR_EXTENSION = 10;
    /**
     * The minimal support for a voted originalCentroid which should occur in at least this many containers, otherwise the
     * extension of the alignment ceases.
     */
    private static final double MINIMAL_SUPPORT = 0.2;
    private int iteration = 0;
    private int numberOfContainers;
    /**
     * Mapping between original container and their already extracted and mapped groups.
     */
    private List<Map<GroupContainer, GroupContainer>> alignedClusters;

    public MultipleStructureAlignment() {
        this.svdSuperimposer = new SVDSuperimposer();
    }

    public void align(List<GroupContainer> containers, int... residueNumbers) {
        this.alignedClusters = new ArrayList<>();
        this.numberOfContainers = containers.size();

        logger.info("creating alignment of {} containers", containers.size());
        logger.info("reference residues {}", Arrays.toString(residueNumbers));

        // move containers to field and extract alignment seed
        Map<GroupContainer, GroupContainer> extractedAlignmentSeeds = containers.parallelStream()
                .collect(Collectors.toMap(Function.identity(), container -> {
                    // select seeds
                    List<Group> extractedGroups = Selection.on(container)
                            .aminoAcids()
                            .residueNumber(residueNumbers)
                            .asFilteredGroups()
                            .collect(Collectors.toList());

                    // remove from parent container
                    container.getGroups().removeAll(extractedGroups);

                    // move to separate container
                    Chain chain = new Chain(extractedGroups);
                    chain.setIdentifier(container.getIdentifier());
                    return chain;
                }));

        GroupContainer consensus = composeConsensus(extractedAlignmentSeeds);

        // align containers with respect to extract groups
        alignWithRespectToConsensus(consensus, extractedAlignmentSeeds);

        alignedClusters.add(extractedAlignmentSeeds);

        extendFragments();
    }

    private GroupContainer composeConsensus(Map<GroupContainer, GroupContainer> containers) {
        // trivial: choose the first entry
//        return containers.entrySet().stream().findFirst().get().getValue();

        logger.debug("composing consensus of {} containers", containers.size());

        // build average coordinate container of everything
        GroupContainer reference = containers.entrySet().stream().findFirst().get().getValue();
        containers.entrySet().stream()
                .map(Map.Entry::getValue)
                .forEach(toMerge -> AbstractConsensusComposer.mergeContainersByCentroid(reference, toMerge));
        return reference;
    }

    public List<GroupContainer> getAlignedContainers() {
        return alignedClusters.stream()
                .map(Map::entrySet)
                .flatMap(Collection::stream)
                .map(Map.Entry::getValue)
                .collect(Collectors.toList());
    }

    private void alignWithRespectToConsensus(GroupContainer consensus, Map<GroupContainer, GroupContainer> containers) {
        containers.entrySet().parallelStream().forEach(entry -> {
            GroupContainer candidate = entry.getValue();
            AlignmentResult alignmentResult = svdSuperimposer.align(consensus, candidate);
            alignmentResult.transform(candidate);
            alignmentResult.transform(entry.getKey());
        });
    }

    public List<Map<GroupContainer, GroupContainer>> getAlignedClusters() {
        return alignedClusters;
    }

    private void extendFragments() {
        iteration++;
        List<Map<GroupContainer, GroupContainer>> newAlignedClusters = new ArrayList<>();
        Map<GroupContainer, GroupContainer> keysToMove = new HashMap<>();
        List<Double> rmsds = new ArrayList<>();

        IntStream.range(0, alignedClusters.size())
                .limit(1)
                .forEach(alignedClusterIndex -> {
            final Map<GroupContainer, GroupContainer> alignedContainers = alignedClusters.get(alignedClusterIndex);
            // determine set of closest groups
            Map<GroupContainer, List<Group>> closestGroups = alignedContainers.entrySet().parallelStream()
                    .collect(Collectors.toMap(Map.Entry::getKey, this::findClosestGroups));

            double[] closestCentroid = LinearAlgebra3D.divide(closestGroups.entrySet().parallelStream()
                    .map(Map.Entry::getValue)
                    .map(list -> list.get(0))
                    .map(LinearAlgebraAtom::centroid)
                    .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), closestGroups.size());

            Map<GroupContainer, Group> selectedGroup = closestGroups.entrySet().parallelStream()
                    .collect(Collectors.toMap(Map.Entry::getKey, entry -> entry.getValue().stream().sorted(Comparator.comparingDouble(group -> squaredDistance(group, closestCentroid))).findFirst().get()));

            selectedGroup.entrySet().forEach(entry -> {
                logger.trace("selected and moving group {}", entry.getValue().getIdentifier());
                GroupContainer origin = entry.getKey();
                Group groupToMove = entry.getValue();
                origin.getGroups().remove(groupToMove);
                alignedContainers.get(origin).getGroups().add(groupToMove);
            });

            GroupContainer consensus = composeConsensus(alignedContainers);

            // align containers with respect to extract groups
            alignedContainers.entrySet().parallelStream().forEach(entry -> {
                        GroupContainer candidate = entry.getValue();
                        AlignmentResult alignmentResult = svdSuperimposer.align(consensus, candidate);
                        alignmentResult.transform(candidate);
                        alignmentResult.transform(entry.getKey());
                        double rmsd = alignmentResult.getAlignmentScore();
                        logger.trace("rmsd of {} and {}: {}", consensus.getIdentifier(),
                                candidate.getIdentifier(),
                                rmsd);

                        // remember freakish occurrences
                        if(rmsd > RMSD_CUTOFF) {
                            keysToMove.put(entry.getKey(), entry.getValue());
                        } else {
                            rmsds.add(rmsd);
                        }
                    });

            keysToMove.forEach(alignedContainers::remove);
            newAlignedClusters.add(alignedContainers);
        });

        // handle freak cluster
        if(alignedClusters.size() > 1) {
            alignedClusters.get(1).entrySet().forEach(entry -> keysToMove.put(entry.getKey(), entry.getValue()));
        }
        GroupContainer freakConsensus = composeConsensus(keysToMove);
        alignWithRespectToConsensus(freakConsensus, keysToMove);
        newAlignedClusters.add(keysToMove);

        alignedClusters = newAlignedClusters;
        int maximalClusterSize = alignedClusters.get(0).size();


        double currentSupport = maximalClusterSize / (double) numberOfContainers;

        logger.debug("iteration: {}, rmsd: {}, maximal support: {}, # clusters: {}",
                iteration,
                rmsds.stream()
                        .mapToDouble(Double::new)
                        .average()
                        .orElseThrow(IllegalArgumentException::new),
                currentSupport,
                alignedClusters.size());

        if(currentSupport > MINIMAL_SUPPORT) {
            extendFragments();
        }
    }

    private List<Group> findClosestGroups(Map.Entry<GroupContainer, GroupContainer> entry) {
        double[] centroid = LinearAlgebra3D.divide(entry.getValue().groups()
                .map(LinearAlgebraAtom::centroid)
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), entry.getValue().getGroups().size());

        return entry.getKey().aminoAcids()
                .sorted(Comparator.comparingDouble(group -> squaredDistance(group, centroid)))
                .limit(NUMBER_OF_GROUPS_AS_CANDIDATES_FOR_EXTENSION)
                .collect(Collectors.toList());
    }

    private double squaredDistance(Group group, double[] centroid) {
        return LinearAlgebra3D.distanceFast(LinearAlgebraAtom.centroid(group), centroid);
    }
}

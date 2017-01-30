package de.bioforscher.jstructure.alignment.multiple;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * Performs a multiple structure alignment.
 * Created by bittrich on 1/25/17.
 */
public class MultipleStructureAlignment {
    private static final Logger logger = LoggerFactory.getLogger(MultipleStructureAlignment.class);
    private final SVDSuperimposer svdSuperimposer;
    private final ClusteringStrategy clusteringStrategy;
    /**
     * The cutoff when a binding site is considered to not be able to extend any further.
     */
    private static final double RMSD_CUTOFF = 1.0;
    private static final double RMSD_CUTOFF_SQUARED = RMSD_CUTOFF * RMSD_CUTOFF;
    private static final int NUMBER_OF_GROUPS_AS_CANDIDATES_FOR_EXTENSION = 10;
    /**
     * The minimal support for a voted originalCentroid which should occur in at least this many containers, otherwise the
     * extension of the alignment ceases.
     */
    private static final double MINIMAL_SUPPORT = 0.8;
    private int numberOfContainers;
    /**
     * Mapping between original container and their already extracted and mapped groups.
     */
    private Map<GroupContainer, GroupContainer> alignedContainers;
    private int maximalAlignmentSize;

    public MultipleStructureAlignment() {
        this.svdSuperimposer = new SVDSuperimposer();
        this.clusteringStrategy = new NaiveClusteringStrategy();
    }

    public List<GroupContainer> align(List<GroupContainer> containers, int... residueNumbers) {
        this.numberOfContainers = containers.size();
        this.maximalAlignmentSize = residueNumbers.length;
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

                    logger.trace("extracted groups for {}: {}", container.getIdentifier(), extractedGroups.stream()
                            .map(Group::getIdentifier)
                            .collect(Collectors.joining(", ", "[", "]")));

                    Chain chain = new Chain(extractedGroups);
                    chain.setIdentifier(container.getIdentifier());
                    return chain;
                }));

        GroupContainer reference = clusteringStrategy.composeConsensus(extractedAlignmentSeeds.entrySet().stream()
                .map(Map.Entry::getValue)
                .collect(Collectors.toList()));

        // align containers with respect to extract groups
        alignedContainers = extractedAlignmentSeeds.entrySet().parallelStream()
                .map(entry -> {
                    GroupContainer candidate = entry.getValue();

                    logger.trace("pairing {} and {}", reference.getGroups().stream()
                            .map(Group::getIdentifier)
                            .collect(Collectors.joining(", ", "[", "]")),
                            candidate.getGroups().stream()
                            .map(Group::getIdentifier)
                            .collect(Collectors.joining(", ", "[", "]")));

                    AlignmentResult alignmentResult = svdSuperimposer.align(reference, candidate);
                    alignmentResult.transform(candidate);
                    alignmentResult.transform(entry.getKey());
                    return entry;
                })
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        extendFragments(1);

        return alignedContainers.entrySet().stream()
                .map(Map.Entry::getValue)
                .collect(Collectors.toList());
    }

    private void extendFragments(int iteration) {
        // determine set of closest groups
        Map<GroupContainer, List<Group>> closestGroups = alignedContainers.entrySet().parallelStream()
                .collect(Collectors.toMap(Map.Entry::getKey, this::findClosestGroups));

        closestGroups.entrySet().forEach(entry -> logger.trace("closest groups for {} are {}",
                    entry.getKey().getIdentifier(),
                    entry.getValue().stream()
                            .map(Group::getIdentifier)
                            .collect(Collectors.joining(", ", "[", "]"))
                ));

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

        GroupContainer reference = clusteringStrategy.composeConsensus(alignedContainers.entrySet().stream()
                .map(Map.Entry::getValue)
                .collect(Collectors.toList()));

        List<Double> rmsds = new ArrayList<>();
        // align containers with respect to extract groups
        alignedContainers = alignedContainers.entrySet().parallelStream()
                .map(entry -> {
                    GroupContainer candidate = entry.getValue();
                    AlignmentResult alignmentResult = svdSuperimposer.align(reference, candidate);
                    alignmentResult.transform(candidate);
                    alignmentResult.transform(entry.getKey());
                    rmsds.add(alignmentResult.getAlignmentScore());
                    logger.trace("rmsd of {} and {}: {}", reference.getIdentifier(),
                            candidate.getIdentifier(),
                            alignmentResult.getAlignmentScore());
                    return entry;
                })
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        logger.debug("average rmsd of iteration {}: {}",
                iteration,
                rmsds.stream().mapToDouble(Double::new).average().getAsDouble());

        if(iteration < 10) {
            extendFragments(iteration + 1);
        }
    }

    public static Collector<GroupContainer, ConsensusComposer, GroupContainer> toConsensus() {
        return Collector.of(ConsensusComposer::new,
                ConsensusComposer::accept,
                ConsensusComposer::combine,
                ConsensusComposer::getConsensus,
                Collector.Characteristics.CONCURRENT);
    }

    static class ConsensusComposer implements Consumer<GroupContainer> {
        @Override
        public void accept(GroupContainer container) {

        }

        ConsensusComposer combine(ConsensusComposer other) {
            return this;
        }

        GroupContainer getConsensus() {
            return null;
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

    /**
     * The specification of strategies to compose or derive one consensus representation of a collection of
     * {@link GroupContainer} objects.
     * Created by S on 28.01.2017.
     */
    interface ClusteringStrategy {
        /**
         * Merges a multitude of structure into one consensus structure which should reasonably well capable of representing
         * this cluster.
         * @param containers the ensemble to create the consensus upon
         * @return a map containing the consensus container as keys and all associated (original) container as value
         */
        Map<GroupContainer, List<GroupContainer>> composeConsensus(List<GroupContainer> containers);
    }

    /**
     * The 'strategy' of creating a consensus by randomly selecting one representative structure.
     * Created by S on 28.01.2017.
     */
    class NaiveClusteringStrategy implements ClusteringStrategy {
        @Override
        public Map<GroupContainer, List<GroupContainer>> composeConsensus(List<GroupContainer> containers) {
            Map<GroupContainer, List<GroupContainer>> map = new HashMap<>();
            map.put(containers.get(0), containers);
            return map;
        }
    }

    class NaiveCentroidClusteringStrategy implements ClusteringStrategy {
        @Override
        public Map<GroupContainer, List<GroupContainer>> composeConsensus(List<GroupContainer> containers) {
            GroupContainer reference = containers.remove(0);
//            svdSuperimposer.align(reference, candidate);
            return null;
        }
    }

    private Atom mergeAtoms(Atom atom1, Atom atom2) {
        atom1.setCoordinates(LinearAlgebra3D.divide(LinearAlgebra3D.add(atom1.getCoordinates(), atom2.getCoordinates()), 2));
        return atom1;
    }

    private GroupContainer mergeContainers(GroupContainer reference, GroupContainer candidate) {
        // get intersecting pairs of atoms wrapped in an atom container
        // 12/14/16 - moved to backboneOnly flag
        Pair<GroupContainer, GroupContainer> containerPair =
                LinearAlgebraAtom.comparableGroupContainerPair(reference,
                        candidate,
                        AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES,
                        AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES);

//        return Combinatorics.sequentialPairsOf(containerPair.getLeft().getAtoms(), containerPair.getRight().getAtoms())
//                .map(AbstractConsensusComposer::mergeAtomPair)
//                .collect(StructureCollectors.toAtomContainer());
        List<Atom> mergedAtoms = new ArrayList<>();
        for(int i = 0; i < containerPair.getLeft().getAtoms().size(); i++) {
            mergedAtoms.add(mergeAtoms(containerPair.getLeft().getAtoms().get(i), containerPair.getRight().getAtoms().get(i)));
        }
        return new Chain(mergedAtoms);
    }
}

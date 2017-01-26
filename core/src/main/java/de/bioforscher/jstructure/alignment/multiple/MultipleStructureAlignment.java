package de.bioforscher.jstructure.alignment.multiple;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
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
//    private Map<String, GroupContainer> alignedFragments;
    /**
     * Mapping between original container and their already extracted and mapped groups.
     */
    private Map<GroupContainer, GroupContainer> alignedContainers;
    private int maximalAlignmentSize;

    public MultipleStructureAlignment() {
        this.svdSuperimposer = new SVDSuperimposer();
    }

    public List<GroupContainer> align(List<GroupContainer> containers, int... residueNumbers) {
        this.numberOfContainers = containers.size();
        this.maximalAlignmentSize = residueNumbers.length;
        logger.info("creating alignment of {} containers", containers.size());
        System.out.println("reference residues " + Arrays.toString(residueNumbers));

        // extract random reference structure
//        GroupContainer referenceStructure = containers.get(0);
//        GroupContainer referenceMotif = Selection.on(referenceStructure)
//                .aminoAcids()
//                .residueNumber(residueNumbers)
//                .cloneElements()
//                .asGroupContainer();

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

        // align all other structure to the reference with respect to the 'seed' of the structure alignment
//        alignedFragments = containers.parallelStream()
//                .map(groupContainer -> {
//                    GroupContainer fragment = Selection.on(groupContainer)
//                            .aminoAcids()
//                            .residueNumber(residueNumbers)
//                            .cloneElements()
//                            .asGroupContainer();
//                    fragment.setIdentifier(groupContainer.getIdentifier());
//                    AlignmentResult alignmentResult = svdSuperimposer.align(referenceMotif, fragment);
//                    alignmentResult.transform(fragment);
//                    alignmentResult.transform(groupContainer);
//                    return fragment;
//                })
//                .collect(Collectors.toMap(GroupContainer::getIdentifier, Function.identity()));

        GroupContainer reference = extractedAlignmentSeeds.entrySet().stream().findFirst().get().getValue();

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

        extendFragments();

        // reassign names to the containers and align them relative to the initial reference motif
//        return alignedFragments.entrySet().stream()
//                .map(Map.Entry::getValue)
//                .collect(Collectors.toList());

        return alignedContainers.entrySet().stream()
                .map(Map.Entry::getValue)
                .collect(Collectors.toList());
    }

    private void extendFragments() {
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
            logger.trace("selected and moving group is {}", entry.getValue().getIdentifier());
            GroupContainer origin = entry.getKey();
            Group groupToMove = entry.getValue();
            origin.getGroups().remove(groupToMove);
            alignedContainers.get(origin).getGroups().add(groupToMove);
        });


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

//    private void extendFragments(List<GroupContainer> candidatesForExtension) {
//        // compute centroid of all previously aligned groups of the candidates
//        //TODO filtering for alignedFragments of maximal size should be faster
//        double[] centroid = LinearAlgebra3D.divide(alignedFragments.entrySet().parallelStream()
////                .filter(entry -> candidatesForExtension.stream().map(GroupContainer::getIdentifier).anyMatch(identifier -> identifier.equals(entry.getKey())))
//                .map(Map.Entry::getValue)
//                .filter(groups -> groups.getGroups().size() == maximalAlignmentSize)
//                .map(LinearAlgebraAtom::centroid)
//                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), candidatesForExtension.size());
////        System.out.println("previous centroid: " + Arrays.toString(centroid));
//
//        List<Group> closestGroups = candidatesForExtension.parallelStream()
//                .map(candidate -> findClosestGroupToCentroid(candidate, centroid))
//                .collect(Collectors.toList());
//        // 'vote' for the next extension of the alignment by majority
//        double[] votedCentroid = LinearAlgebra3D.divide(closestGroups.stream()
//                .map(LinearAlgebraAtom::centroid)
//                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), closestGroups.size());
////        System.out.println("voted centroid: " + Arrays.toString(votedCentroid));
//
//        // for each candidate: find group closest to voted centroid
//        Map<GroupContainer, Group> selectedCandidatesForExtension = new HashMap<>();
//        candidatesForExtension.parallelStream().forEach(candidate -> {
//                    Group closestGroup = findClosestGroupToCentroid(candidate, votedCentroid);
//                    double squaredDistanceToCentroid = LinearAlgebra3D.distanceFast(LinearAlgebraAtom.centroid(closestGroup), votedCentroid);
//                    if(squaredDistanceToCentroid < RMSD_CUTOFF_SQUARED) {
//                        selectedCandidatesForExtension.put(candidate, closestGroup);
//                    }
//                });
//
//        if(selectedCandidatesForExtension.size() == 0 || selectedCandidatesForExtension.size() < numberOfContainers * MINIMAL_SUPPORT) {
//            System.out.println(selectedCandidatesForExtension.size());
//            System.out.println(numberOfContainers);
//            System.out.println(numberOfContainers * MINIMAL_SUPPORT);
//            return;
//        }
//
//        // extend fragments
//        selectedCandidatesForExtension.entrySet().parallelStream().forEach(candidateEntry -> {
//            alignedFragments.get(candidateEntry.getKey().getIdentifier()).getGroups().add(candidateEntry.getValue());
//        });
//        // select new reference motif
//        GroupContainer referenceMotif = alignedFragments.get(selectedCandidatesForExtension.entrySet().stream().findFirst().get().getKey().getIdentifier());
//        // update orientation
//
//
//        // start next iteration with extended candidates
//        extendFragments(selectedCandidatesForExtension.entrySet().stream()
//                .map(Map.Entry::getKey)
//                .collect(Collectors.toList()));
//    }
//
//    private Group findClosestGroupToCentroid(GroupContainer candidate, double[] extensionCentroid) {
//        List<Group> handledGroupsForThisContainer = alignedFragments.get(candidate.getIdentifier()).getGroups();
//        return candidate.aminoAcids()
//                .filter(group -> !handledGroupsForThisContainer.contains(group))
//                .min((g1, g2) -> {
//                    double[] centroid1 = LinearAlgebraAtom.centroid(g1);
//                    double[] centroid2 = LinearAlgebraAtom.centroid(g2);
//                    double distance1 = LinearAlgebra3D.distanceFast(centroid1, extensionCentroid);
//                    double distance2 = LinearAlgebra3D.distanceFast(centroid2, extensionCentroid);
//                    return Double.compare(distance1, distance2);
//                })
//                .orElseThrow(() -> new IllegalArgumentException("candidate is empty: " + candidate + System.lineSeparator() + candidate.composePDBRecord()));
//    }
}

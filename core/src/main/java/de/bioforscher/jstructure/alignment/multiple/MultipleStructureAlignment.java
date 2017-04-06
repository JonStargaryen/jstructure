package de.bioforscher.jstructure.alignment.multiple;

import de.bioforscher.jstructure.alignment.StructureAlignmentResult;
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

/**
 * Performs a multiple structure alignment.
 * TODO parallelisms
 * Created by bittrich on 1/25/17.
 */
@Deprecated
public class MultipleStructureAlignment {
    private static final Logger logger = LoggerFactory.getLogger(MultipleStructureAlignment.class);
    private final SVDSuperimposer svdSuperimposer;
    /**
     * The cutoff when a binding site is considered to not be able to extend any further.
     */
    private static final double RMSD_CUTOFF = 4.0;
    private static final int NUMBER_OF_GROUPS_AS_CANDIDATES_FOR_EXTENSION = 10;
    /**
     * The minimal support for a voted originalCentroid which should occur in at least this many containers, otherwise the
     * extension of the alignment ceases.
     */
    private static final double MINIMAL_SUPPORT = 0.8;
    private int iteration = 0;
    private int numberOfContainers;
    /*
     * Mapping between original container and their already extracted and mapped groups. Key: original container, value:
     * extracted fragments representing the original container in the alignment. Groups present in the latter will be
     * removed from the first.
     */
    /**
     * The collection of containers which are still extending.
     */
    private Map<GroupContainer, GroupContainer> containersExtending;
    /**
     * The collection of containers which cannot be extended any further as their observed RMSD was above the threshold.
     */
    private Map<GroupContainer, GroupContainer> containersFinished;

    public MultipleStructureAlignment() {
        svdSuperimposer = new SVDSuperimposer();
        containersFinished = new HashMap<>();
    }

    public void align(List<GroupContainer> containers, int... residueNumbers) {
        numberOfContainers = containers.size();

        logger.info("creating alignment of {} containers", containers.size());
        logger.info("reference residues {}", Arrays.toString(residueNumbers));

        // move containers to field and extract alignment seed
        containersExtending = containers.stream()
                .collect(Collectors.toMap(Function.identity(), container -> {
                    // select seeds
                    List<Group> extractedGroups = Selection.on(container)
                            .aminoAcids()
                            .residueNumber(residueNumbers)
                            .cloneElements()
                            .asFilteredGroups()
                            .collect(Collectors.toList());

                    // remove from parent container
                    container.getGroups().removeAll(extractedGroups);

                    // move to separate container - 'remember name'
                    Chain chain = new Chain(extractedGroups);
                    chain.setIdentifier(container.getIdentifier());
                    return chain;
                }));

        // align containers with respect to extract groups
        GroupContainer consensus = composeConsensus(containersExtending);
        alignWithRespectToConsensus(consensus, containersExtending);

        // start iterative loop
        extendFragments();

        // move previously extending containers to finished category
        containersExtending.entrySet().forEach(entry -> containersFinished.put(entry.getKey(), entry.getValue()));

        // align everything with respect to the smallest container
        int minimalSize = containersFinished.entrySet().stream()
                .map(Map.Entry::getValue)
                .map(GroupContainer::getGroups)
                .mapToInt(Collection::size)
                .min()
                .orElseThrow(() -> new IllegalArgumentException("container list is empty"));

        // extract containers of minimal size and compose consensus of them
        Map<GroupContainer, GroupContainer> minimalSizedContainers = containersFinished.entrySet().stream()
                .filter(entry -> entry.getValue().getGroups().size() == minimalSize)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        GroupContainer minimalConsensus = composeConsensus(minimalSizedContainers);

        // align everything with respect to consensus
        alignWithRespectToConsensus(minimalConsensus, containersFinished);

        // reassign extracted residues to original containers
        containersFinished.entrySet().parallelStream().forEach(entry -> {
            GroupContainer originalContainer = entry.getKey();
            GroupContainer extractedGroups = entry.getValue();
            originalContainer.getGroups().addAll(extractedGroups.getGroups());
        });
    }

    private GroupContainer composeConsensus(Map<GroupContainer, GroupContainer> containers) {
        //TODO maybe there is a more robust/sophisticated/meaningful/reasonable consensus building strategy
        // build average coordinate container of everything
        GroupContainer reference = containers.entrySet().stream()
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException("collection to create consensus upon is empty"))
                .getValue();
        containers.entrySet().stream()
                .map(Map.Entry::getValue)
                .forEach(toMerge -> AbstractConsensusComposer.mergeContainersByCentroid(reference, toMerge));
        return reference;
    }

    private void alignWithRespectToConsensus(GroupContainer consensus, Map<GroupContainer, GroupContainer> containers) {
        containers.entrySet().parallelStream().forEach(entry -> {
            GroupContainer originalGroups = entry.getKey();
            GroupContainer extractedGroups = entry.getValue();
            // align extracted groups to consensus
            StructureAlignmentResult alignmentResult = svdSuperimposer.align(consensus, extractedGroups);

            // employ alignment - TODO this could more efficiently happen in-place
            alignmentResult.transform(originalGroups);
            alignmentResult.transform(extractedGroups);
        });
    }

    private void extendFragments() {
        iteration++;

        // determine set of closest groups
        Map<GroupContainer, List<Group>> closestGroups = containersExtending.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, this::findClosestGroups));

        double[] closestCentroid = LinearAlgebra3D.divide(closestGroups.entrySet().stream()
                .map(Map.Entry::getValue)
                .map(list -> list.get(0))
                .map(LinearAlgebraAtom::centroid)
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), closestGroups.size());

        // select group closest to centroid
        Map<GroupContainer, Group> selectedGroup = closestGroups.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> entry.getValue().stream()
                        .sorted(Comparator.comparingDouble(group -> calculateSquaredDistance(group, closestCentroid)))
                        .findFirst()
//                        .orElseThrow(() -> new IllegalArgumentException("did not find nearby group"))));
                        .get()));

        selectedGroup.entrySet().forEach(entry -> {
            GroupContainer originalGroups = entry.getKey();
            GroupContainer extractedGroups = containersExtending.get(originalGroups);
            Group extractedGroup = entry.getValue();

            // move selected group from original to extracted groups
            logger.trace("selected and moving group {}", extractedGroup.getIdentifier());
            // fix 2/2/17 - .getGroups() points to
//            originalGroups.getGroups().remove(extractedGroup);
            //TODO generic removeFromParent-function? - this is really fragile and can break easily
            extractedGroup.getParentChain().getGroups().remove(extractedGroup);
            extractedGroups.getGroups().add(extractedGroup);
        });

        GroupContainer consensus = composeConsensus(containersExtending);
        logger.trace("current consensus\n{}", consensus.composePDBRecord());
        List<GroupContainer> containersCeasedExtending = new ArrayList<>();

        // align containers with respect to extract groups
        containersExtending.entrySet().parallelStream().forEach(entry -> {
            GroupContainer originalGroups = entry.getKey();
            GroupContainer extractedGroups = entry.getValue();
            StructureAlignmentResult alignmentResult = svdSuperimposer.align(consensus, extractedGroups);
            alignmentResult.transform(originalGroups);
            alignmentResult.transform(extractedGroups);
            double rmsd = alignmentResult.getAlignmentScore();
            logger.trace("rmsd of {} and {}: {}", consensus.getIdentifier(),
                    extractedGroups.getIdentifier(),
                    rmsd);

            // container is dissimilar: track as 'finished', i.e. it can not be extended any further
            if(rmsd > RMSD_CUTOFF) {
                containersCeasedExtending.add(originalGroups);
                containersFinished.put(originalGroups, extractedGroups);
            }
        });

        // remove containers ceased extending from list of still considered candidates for extension
        containersCeasedExtending.forEach(containersExtending::remove);

        int maximalClusterSize = containersExtending.size();
        double currentSupport = maximalClusterSize / (double) numberOfContainers;

        logger.debug("iteration: {}, current support: {}, containers: {} + {} = {}",
                iteration,
                currentSupport,
                containersExtending.size(),
                containersFinished.size(),
                containersExtending.size() + containersFinished.size());

        if(currentSupport > MINIMAL_SUPPORT) {
            extendFragments();
        }
    }

    private List<Group> findClosestGroups(Map.Entry<GroupContainer, GroupContainer> entry) {
        double[] centroid = LinearAlgebra3D.divide(entry.getValue().groups()
                .map(LinearAlgebraAtom::centroid)
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), entry.getValue().getGroups().size());

        return entry.getKey().aminoAcids()
                .sorted(Comparator.comparingDouble(group -> calculateSquaredDistance(group, centroid)))
                .limit(NUMBER_OF_GROUPS_AS_CANDIDATES_FOR_EXTENSION)
                .collect(Collectors.toList());
    }

    private double calculateSquaredDistance(Group group, double[] centroid) {
        return LinearAlgebra3D.distanceFast(LinearAlgebraAtom.centroid(group), centroid);
    }

    public Map<GroupContainer, GroupContainer> getAlignedContainerMap() {
        return containersFinished;
    }

    public List<GroupContainer> getAlignedFragments() {
        return containersFinished.entrySet().stream()
                .map(Map.Entry::getValue)
                .collect(Collectors.toList());
    }
}

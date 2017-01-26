package de.bioforscher.jstructure.alignment.bindingsite;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 1/25/17.
 */
public class MultipleStructureAlignment {
    private final SVDSuperimposer svdSuperimposer;
    /**
     * The cutoff when a binding site is considered to not be able to extend any further.
     */
    private static final double RMSD_CUTOFF = 1.0;
    private static final double RMSD_CUTOFF_SQUARED = RMSD_CUTOFF * RMSD_CUTOFF;
    /**
     * The minimal support for a voted originalCentroid which should occur in at least this many containers, otherwise the
     * extension of the alignment ceases.
     */
    private static final double MINIMAL_SUPPORT = 0.8;
    private int numberOfContainers;
    private Map<String, List<Group>> alignedFragments;
    private int maximalAlignmentSize;

    public MultipleStructureAlignment() {
        this.svdSuperimposer = new SVDSuperimposer();
    }

    public List<GroupContainer> align(List<GroupContainer> containers, int... residueNumbers) {
        this.numberOfContainers = containers.size();
        this.maximalAlignmentSize = residueNumbers.length;

        // extract random reference structure
        GroupContainer referenceStructure = containers.get(0);
        GroupContainer referenceMotif = Selection.on(referenceStructure)
                .aminoAcids()
                .residueNumber(residueNumbers)
                .cloneElements()
                .asGroupContainer();

        // align all other structure to the reference with respect to the 'seed' of the structure alignment
        alignedFragments = containers.parallelStream()
                .map(groupContainer -> {
                    GroupContainer fragment = Selection.on(groupContainer)
                            .aminoAcids()
                            .residueNumber(residueNumbers)
                            .cloneElements()
                            .asGroupContainer();
                    fragment.setIdentifier(groupContainer.getIdentifier());
                    AlignmentResult alignmentResult = svdSuperimposer.align(referenceMotif, fragment);
                    alignmentResult.transform(fragment);
                    alignmentResult.transform(groupContainer);
                    return fragment;
                })
                .collect(Collectors.toMap(GroupContainer::getIdentifier, GroupContainer::getGroups));

        extendFragments(containers);

        // reassign names to the containers and align them relative to the initial reference motif
        return alignedFragments.entrySet().stream()
                .map(entry -> {
                    Chain container = new Chain(entry.getValue());
                    container.setIdentifier(entry.getKey());
                    //TODO needed?
//                    svdSuperimposer.align(referenceMotif, container).transform(container);
                    return container;
                })
                .collect(Collectors.toList());
    }

    private void extendFragments(List<GroupContainer> candidatesForExtension) {
        // compute centroid of all previously aligned groups of the candidates
        //TODO filtering for alignedFragments of maximal size should be faster
        double[] centroid = LinearAlgebra3D.divide(alignedFragments.entrySet().parallelStream()
//                .filter(entry -> candidatesForExtension.stream().map(GroupContainer::getIdentifier).anyMatch(identifier -> identifier.equals(entry.getKey())))
                .map(Map.Entry::getValue)
                .filter(groups -> groups.size() == maximalAlignmentSize)
                .map(Chain::new)
                .map(LinearAlgebraAtom::centroid)
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), candidatesForExtension.size());
//        System.out.println("previous centroid: " + Arrays.toString(centroid));

        List<Group> closestGroups = candidatesForExtension.parallelStream()
                .map(candidate -> findClosestGroupToCentroid(candidate, centroid))
                .collect(Collectors.toList());
        // 'vote' for the next extension of the alignment by majority
        double[] votedCentroid = LinearAlgebra3D.divide(closestGroups.stream()
                .map(LinearAlgebraAtom::centroid)
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), closestGroups.size());
//        System.out.println("voted centroid: " + Arrays.toString(votedCentroid));

        // for each candidate: find group closest to voted centroid
        Map<GroupContainer, Group> selectedCandidatesForExtension = new HashMap<>();
        candidatesForExtension.parallelStream().forEach(candidate -> {
                    Group closestGroup = findClosestGroupToCentroid(candidate, votedCentroid);
                    double squaredDistanceToCentroid = LinearAlgebra3D.distanceFast(LinearAlgebraAtom.centroid(closestGroup), votedCentroid);
                    if(squaredDistanceToCentroid < RMSD_CUTOFF_SQUARED) {
                        selectedCandidatesForExtension.put(candidate, closestGroup);
                    }
                });

        if(selectedCandidatesForExtension.size() == 0 || selectedCandidatesForExtension.size() < numberOfContainers * MINIMAL_SUPPORT) {
            System.out.println(selectedCandidatesForExtension.size());
            System.out.println(numberOfContainers);
            System.out.println(numberOfContainers * MINIMAL_SUPPORT);
            return;
        }

        // extend fragments
        selectedCandidatesForExtension.entrySet().parallelStream().forEach(candidateEntry -> {
            alignedFragments.get(candidateEntry.getKey().getIdentifier()).add(candidateEntry.getValue());
        });
        // select new reference motif
        GroupContainer referenceMotif = new Chain(alignedFragments.get(selectedCandidatesForExtension.entrySet().stream().findFirst().get().getKey().getIdentifier()));
        // update orientation


        // start next iteration with extended candidates
        extendFragments(selectedCandidatesForExtension.entrySet().stream()
                .map(Map.Entry::getKey)
                .collect(Collectors.toList()));
    }

    private Group findClosestGroupToCentroid(GroupContainer candidate, double[] extensionCentroid) {
        List<Group> handledGroupsForThisContainer = alignedFragments.get(candidate.getIdentifier());
        return candidate.aminoAcids()
                .filter(group -> !handledGroupsForThisContainer.contains(group))
                .min((g1, g2) -> {
                    double[] centroid1 = LinearAlgebraAtom.centroid(g1);
                    double[] centroid2 = LinearAlgebraAtom.centroid(g2);
                    double distance1 = LinearAlgebra3D.distanceFast(centroid1, extensionCentroid);
                    double distance2 = LinearAlgebra3D.distanceFast(centroid2, extensionCentroid);
                    return Double.compare(distance1, distance2);
                })
                .orElseThrow(() -> new IllegalArgumentException("candidate is empty: " + candidate + System.lineSeparator() + candidate.composePDBRecord()));
    }
}

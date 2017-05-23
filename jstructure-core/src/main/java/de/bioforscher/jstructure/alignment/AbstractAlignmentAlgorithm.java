package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm<R extends AlignmentResult> implements AlignmentAlgorithm<R> {
    private static final Logger logger = LoggerFactory.getLogger(AbstractAlignmentAlgorithm.class);
    protected final Set<String> minimalSetOfAtomNames;
    protected final Set<String> maximalSetOfAtomNames;

    public AbstractAlignmentAlgorithm() {
        this(Collections.emptySet(), AminoAcidFamily.ATOM_NAMES.ALL_ATOM_NAMES);
    }

    public AbstractAlignmentAlgorithm(Set<String> minimalSetOfAtomNames, Set<String> maximalSetOfAtomNames) {
        this.minimalSetOfAtomNames = minimalSetOfAtomNames;
        this.maximalSetOfAtomNames = maximalSetOfAtomNames;
    }

    /**
     * Computes the root-mean-square deviation between 2 atom sets.
     * @param atomContainer1 a collection of atoms
     * @param atomContainer2 another collection of atoms
     * @return the RMSD value of this alignment
     */
    static double calculateRmsd(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        Pair<GroupContainer, GroupContainer> groupContainerPair = comparableGroupContainerPair(atomContainer1, atomContainer2);
        double msd = Combinatorics.sequentialPairsOf(groupContainerPair.getLeft().getAtoms(), groupContainerPair.getRight().getAtoms())
                .mapToDouble(pair -> LinearAlgebra.on(pair.getLeft().getCoordinates()).distanceFast(pair.getRight().getCoordinates()))
                .average()
                .orElseThrow(() -> new IllegalArgumentException("cannot compute rmsd for empty or non-intersecting containers"));
        return Math.sqrt(msd);
    }

    /**
     * Creates a comparable entity of 2 {@link AtomContainer} objects. They fulfill certain criteria:
     * <ul>
     *     <li>{@link Atom}s are cloned
     *     <li>for each {@link Group} only shared atoms are retained</li>
     *     <li>{@link Atom}s are in identical ordering</li>
     *     <li>{@link Atom}s are wrapped in {@link Group} objects if they were initially</li>
     * </ul>
     * @param container1 a collection of reference atoms
     * @param container2 a collection of candidate atoms
     * @param minimalSetOfAtomNames the lower bound of required atom names present (e.g. by setting CA, every time an
     *                              amino acid missing the alpha carbon will result in an exception thrown)
     * @param maximalSetOfAtomNames the upper bound of required atom names present (e.g. by setting CA, everything else
     *                              will be dropped, even when both amino acids would share more atoms)
     * @return a pair of both collections which can now be aligned
     */
    static Pair<GroupContainer, GroupContainer> comparableGroupContainerPair(AtomContainer container1,
                                                                                    AtomContainer container2,
                                                                                    Set<String> minimalSetOfAtomNames,
                                                                                    Set<String> maximalSetOfAtomNames) {
        //TODO what happens for alternative positions?
        GroupContainer groupContainer1 = cloneIntoGroupContainer(container1);
        GroupContainer groupContainer2 = cloneIntoGroupContainer(container2);

        int limitingSize = Math.min(groupContainer1.getGroups().size(), groupContainer2.getGroups().size());

        List<Group> groups1 = new ArrayList<>();
        List<Group> groups2 = new ArrayList<>();
        for(int groupIndex = 0; groupIndex < limitingSize; groupIndex++) {
            Group group1 = groupContainer1.getGroups().get(groupIndex);
            Group group2 = groupContainer2.getGroups().get(groupIndex);
            Pair<List<Atom>, List<Atom>> sharedAtoms = selectSharedAtoms(group1,
                    group2,
                    minimalSetOfAtomNames,
                    maximalSetOfAtomNames);

            // remove additional/not-shared atoms
            group1.getAtoms().clear();
            group1.getAtoms().addAll(sharedAtoms.getLeft());
            group2.getAtoms().clear();
            group2.getAtoms().addAll(sharedAtoms.getRight());

            logger.trace("shared atoms between {} and {}: {}", group1, group2, sharedAtoms);

            groups1.add(group1);
            groups2.add(group2);
        }

        return new Pair<>(new Chain(groups1), new Chain(groups2));
    }

    /**
     * Returns the set of atoms shared by both containers.
     */
    static Pair<GroupContainer, GroupContainer> comparableGroupContainerPair(AtomContainer atomContainer1,
                                                                                    AtomContainer atomContainer2) {
        return comparableGroupContainerPair(atomContainer1,
                atomContainer2,
                Collections.emptySet(),
                Collections.emptySet());
    }

    private static Pair<List<Atom>, List<Atom>> selectSharedAtoms(Group group1,
                                                                  Group group2,
                                                                  Set<String> minimalSetOfAtomNames,
                                                                  Set<String> maximalSetOfAtomNames) {
        // no ordering is enforced - rather the filtering by atom names ensures the identical occurrence of atoms
        Set<String> sharedAtomNames = determineSharedAtomNames(group1,
                group2,
                minimalSetOfAtomNames,
                maximalSetOfAtomNames);

        return new Pair<>(selectAtoms(group1, sharedAtomNames), selectAtoms(group2, sharedAtomNames));
    }

    private static List<Atom> selectAtoms(Group group, Set<String> atomNamesToSelect) {
        // fix 02/08/17 - mere retaining does not ensure correct ordering
        return atomNamesToSelect.stream()
                .flatMap(name -> group.atoms()
                        .filter(atom -> atom.getName().equals(name)))
                .collect(Collectors.toList());
    }

    private static Set<String> determineSharedAtomNames(Group group1,
                                                        Group group2,
                                                        Set<String> minimalSetOfAtomNames,
                                                        Set<String> maximalSetOfAtomNames) {
        // present atom names in each group, determine the shared ones
        Set<String> sharedAtomNames = group1.atoms()
                .map(Atom::getName)
                .collect(Collectors.toSet());
        sharedAtomNames.retainAll(group2.atoms()
                .map(Atom::getName)
                .collect(Collectors.toSet()));

        // fail if the minimal set of atoms is not fulfilled
        if(!sharedAtomNames.containsAll(minimalSetOfAtomNames)) {
            throw new IllegalArgumentException("alignment could not fulfill minimal required atom names" +
                    System.lineSeparator() + group1.getIdentifier() + " atoms: " + group1.atoms().map(Atom::getName).collect(Collectors.joining(", ")) +
                    System.lineSeparator() + group2.getIdentifier() + " atoms: " + group2.atoms().map(Atom::getName).collect(Collectors.joining(", ")) +
                    System.lineSeparator() + "shared: " + sharedAtomNames.stream().collect(Collectors.joining(", ")) +
                    System.lineSeparator() + "required: " + minimalSetOfAtomNames.stream().collect(Collectors.joining(", ")));
        }

        // ignore all atoms not explicitly wanted in the alignment
        if(!maximalSetOfAtomNames.isEmpty()) {
            sharedAtomNames.retainAll(maximalSetOfAtomNames);
        }

        return sharedAtomNames;
    }

    /**
     * Maps an {@link AtomContainer} to a distinct copy of itself.
     * @param atomContainer the container to clone
     * @return a {@link GroupContainer} (i.e. motivated by interfaces handling {@link AtomContainer}) which can readily
     * manipulated without affecting the original entry
     */
    private static GroupContainer cloneIntoGroupContainer(AtomContainer atomContainer) {
        return atomContainer.atoms()
                // map to group level - will still collect 'dangling' atoms into a group
                .map(Atom::getParentGroup)
                .distinct()
                // clone
                .map(Group::new)
                .collect(StructureCollectors.toGroupContainer());
    }
}

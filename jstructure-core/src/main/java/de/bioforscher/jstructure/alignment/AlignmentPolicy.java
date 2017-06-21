package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Rules for alignments.
 * Created by bittrich on 6/19/17.
 */
public interface AlignmentPolicy {
    /**
     *
     */
    @FunctionalInterface
    interface AtomMapping {
        /**
         * Traverses both containers and returns a comparable arrangement of atoms.
         * @param reference the reference container
         * @param query the query container
         * @return a collection of compatible, paired atoms of both containers
         */
        List<Pair<Atom, Atom>> determineAtomMapping(GroupContainer reference, GroupContainer query);
    }

    /**
     * Describes which atoms or groups are used (and in which way) to align both containers.
     * TODO split into validation, selection/filtering and matching
     */
    enum MatchingBehavior implements AtomMapping {
        /**
         * Assumes equal number of atoms.
         */
        ATOM_INDEX((reference, query) -> {
            ensureMatchingAtomCount(reference, query);
            return IntStream.range(0, reference.getAtoms().size())
                    .mapToObj(index -> new Pair<>(reference.getAtoms().get(index), query.getAtoms().get(index)))
                    .collect(Collectors.toList());
        }),
        /**
         * Assumes equal number of groups. Allows for variable matched groups (e.g. Ala vs Ile), they must share at
         * least 1 atom however.
         */
        COMPARABLE_ATOM_NAMES((reference, query) -> {
            ensureMatchingGroupCount(reference, query);
            return IntStream.range(0, reference.getGroups().size())
                    .mapToObj(index -> determineSharedAtoms(reference.getGroups().get(index), query.getGroups().get(index)))
                    .flatMap(Collection::stream)
                    .collect(Collectors.toList());
        }),
        /**
         * Matches amino acids of both containers. Assumes equal number of amino acids. Allows for variable matched
         * groups (e.g. Ala vs Ile), they must share at least 1 atom however.
         */
        AMINO_ACIDS_COMPARABLE_ATOM_NAMES((reference, query) -> {
            GroupContainer referenceAminoAcids = reference.aminoAcids().collect(StructureCollectors.toGroupContainer());
            GroupContainer queryAminoAcids = query.aminoAcids().collect(StructureCollectors.toGroupContainer());
            ensureMatchingGroupCount(referenceAminoAcids, queryAminoAcids);
            return IntStream.range(0, referenceAminoAcids.getGroups().size())
                    .mapToObj(index -> determineSharedAtoms(referenceAminoAcids.getGroups().get(index), queryAminoAcids.getGroups().get(index)))
                    .flatMap(Collection::stream)
//                    .peek(System.out::println)
                    .collect(Collectors.toList());
        }),
        AMINO_ACIDS_COMPARABLE_BACKBONE_ATOM_NAMES((reference, query) -> {
            GroupContainer referenceAminoAcids = reference.aminoAcids().collect(StructureCollectors.toGroupContainer());
            GroupContainer queryAminoAcids = query.aminoAcids().collect(StructureCollectors.toGroupContainer());
            ensureMatchingGroupCount(referenceAminoAcids, queryAminoAcids);
            return IntStream.range(0, referenceAminoAcids.getGroups().size())
                    .mapToObj(index -> determineSharedBackboneAtoms(referenceAminoAcids.getGroups().get(index), queryAminoAcids.getGroups().get(index)))
                    .flatMap(Collection::stream)
//                    .peek(System.out::println)
                    .collect(Collectors.toList());
        }),
        AMINO_ACIDS_ALPHA_CARBONS_TOLERANT((reference, query) -> {
            List<AminoAcid> referenceAminoAcids = reference.aminoAcids().collect(Collectors.toList());
            List<AminoAcid> queryAminoAcids = query.aminoAcids().collect(Collectors.toList());
            int limitingSize = Math.min(referenceAminoAcids.size(), queryAminoAcids.size());
            return IntStream.range(0, limitingSize)
                    //TODO fallback to centroid?
                    .mapToObj(index -> new Pair<>(referenceAminoAcids.get(index).getCa(), queryAminoAcids.get(index).getCa()))
                    .collect(Collectors.toList());
        });

        private final AtomMapping atomMapping;

        MatchingBehavior(AtomMapping atomMapping) {
            this.atomMapping = atomMapping;
        }

        /**
         * Pairs two {@link de.bioforscher.jstructure.model.structure.container.AtomContainer} entities in a comparable
         * way. Will determine the set of shared atom names and pair matching names of both containers.
         * @param referenceGroup the reference container
         * @param queryGroup the query container
         * @return a collection of compatible atom pairs
         */
        private static List<Pair<Atom, Atom>> determineSharedAtoms(Group referenceGroup, Group queryGroup) {
            // determine set of shared atom names
            Set<String> sharedAtomNames = referenceGroup.atoms()
                    .map(Atom::getName)
                    .filter(atomName -> queryGroup.atoms().anyMatch(queryGroupAtom -> atomName.equals(queryGroupAtom.getName())))
                    .collect(Collectors.toSet());

            // validate that both groups share atoms
            if(sharedAtomNames.isEmpty()) {
                throw new IllegalArgumentException("groups " + referenceGroup + " and " + queryGroup + " do not share " +
                        "any atoms at all - cannot align");
            }

            // pair atoms
            return sharedAtomNames.stream()
                    .map(atomName -> new Pair<>(selectAtom(referenceGroup, atomName), selectAtom(queryGroup, atomName)))
                    .collect(Collectors.toList());
        }

        private static List<Pair<Atom, Atom>> determineSharedBackboneAtoms(Group referenceGroup, Group queryGroup) {
            return determineSharedAtoms(referenceGroup, queryGroup).stream()
                    .filter(pair -> AminoAcid.isBackboneAtom(pair.getLeft()) && AminoAcid.isBackboneAtom(pair.getRight()))
                    .collect(Collectors.toList());
        }

        /**
         * Native implementation of selection method. Could use
         * {@link de.bioforscher.jstructure.model.structure.selection.Selection} as well, but performance ought to be
         * favorable for the straight-forward implementation.
         * @param group the container to process
         * @param atomName the atom name to retrieve
         * @return the desired atom
         */
        private static Atom selectAtom(Group group, String atomName) {
            return group.atoms()
                    .filter(atom -> atomName.equals(atom.getName()))
                    //TODO test case for multiple atom names within the same group
                    .findFirst()
                    // presence is guaranteed by looking at shared atom names
                    .get();
        }

        /**
         * Checks for an agreement in total atom count.
         * @param reference the reference container
         * @param query the query container
         * @throws AlignmentException when both containers do not match in size
         */
        private static void ensureMatchingAtomCount(GroupContainer reference, GroupContainer query) {
            if(reference.getAtoms().size() != query.getAtoms().size()) {
                throw new AlignmentException("atom count in both containers does not match! " +
                        reference.getAtoms().size() + " vs " + query.getAtoms().size());
            }
        }

        /**
         * Checks for an agreement in total group count.
         * @param reference the reference container
         * @param query the query container
         * @throws AlignmentException when both containers do not match in size
         */
        private static void ensureMatchingGroupCount(GroupContainer reference, GroupContainer query) {
            if(reference.getGroups().size() != query.getGroups().size()) {
                throw new AlignmentException("group count in both containers does not match! " + System.lineSeparator()
                        + "length: " + reference.getGroups().size() + " vs " + query.getGroups().size() +
                        System.lineSeparator() + "sequences:" + System.lineSeparator() +
                        reference.getAminoAcidSequence() + System.lineSeparator() + query.getAminoAcidSequence());
            }
        }

        @Override
        public List<Pair<Atom, Atom>> determineAtomMapping(GroupContainer reference, GroupContainer query) {
            return atomMapping.determineAtomMapping(reference, query);
        }
    }

    /**
     * Determines whether the original query will be manipulated or if a copy is created whose coordinates will be
     * manipulated.
     */
    enum ManipulationBehavior {
        /**
         * Causes the original reference and query to be manipulated (i.e. centered and the query is rotated onto the
         * reference).
         */
        INPLACE,
        /**
         * Creates copies of the provided reference and query and manipulates those entities. The original arguments will
         * not be modified by any means.
         */
        COPY
    }
}

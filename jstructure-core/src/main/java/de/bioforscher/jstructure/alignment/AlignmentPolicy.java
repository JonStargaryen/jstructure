package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
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
         *
         * @param reference
         * @param query
         * @return
         */
        List<Pair<Atom, Atom>> determineAtomMapping(GroupContainer reference, GroupContainer query);
    }

    /**
     * Describes which atoms or groups are used (and in which way) to align both containers.
     */
    enum MatchingBehavior implements AtomMapping {
        /**
         *
         */
        BY_ATOM_INDEX((r, q) -> {
            ensureMatchingAtomCount(r, q);
            return IntStream.range(0, r.getAtoms().size())
                    .mapToObj(index -> new Pair<>(r.getAtoms().get(index), q.getAtoms().get(index)))
                    .collect(Collectors.toList());
        }),
        BY_COMPARABLE_ATOM_NAMES((r, q) -> {
            ensureMatchingGroupCount(r, q);
            return IntStream.range(0, r.getGroups().size())
                    .mapToObj(index -> determineSharedAtoms(r.getGroups().get(index), q.getGroups().get(index)))
                    .flatMap(Collection::stream)
                    .collect(Collectors.toList());
        });

        private final AtomMapping atomMapping;

        MatchingBehavior(AtomMapping atomMapping) {
            this.atomMapping = atomMapping;
        }

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

        /**
         * Native implementation of selection method. Could use
         * {@link de.bioforscher.jstructure.model.structure.selection.Selection} as well, but performance ought to be
         * favorable for the straight-forward implementation.
         * @param group
         * @param atomName
         * @return
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
         *
         * @param reference
         * @param query
         */
        private static void ensureMatchingAtomCount(GroupContainer reference, GroupContainer query) {
            if(reference.getAtoms().size() != query.getAtoms().size()) {
                throw new AlignmentException("atom count in both containers does not match! " +
                        reference.getAtoms().size() + " vs " + query.getAtoms().size());
            }
        }

        /**
         *
         * @param reference
         * @param query
         */
        private static void ensureMatchingGroupCount(GroupContainer reference, GroupContainer query) {
            if(reference.getGroups().size() != query.getGroups().size()) {
                throw new AlignmentException("group count in both containers does not match! " +
                        reference.getGroups().size() + " vs " + query.getGroups().size());
            }
        }

        @Override
        public List<Pair<Atom, Atom>> determineAtomMapping(GroupContainer reference, GroupContainer query) {
            return atomMapping.determineAtomMapping(reference, query);
        }
    }

    /**
     *
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

package de.bioforscher.jstructure.align;

import de.bioforscher.jstructure.model.structure.container.GroupContainer;

/**
 * Builder for rather complex structure alignments.
 * Created by S on 12.07.2017.
 */
public class StructureAlignmentQuery {
    private final GroupContainer reference;
    private final GroupContainer query;
    private final AlignmentPolicy.AtomMapping atomMapping;

    private StructureAlignmentQuery(MatchingBehaviorStep matchingBehaviorStep) {
        this.reference = matchingBehaviorStep.reference;
        this.query = matchingBehaviorStep.query;
        this.atomMapping = matchingBehaviorStep.atomMapping;
    }

    public GroupContainer getReference() {
        return reference;
    }

    public GroupContainer getQuery() {
        return query;
    }

    public AlignmentPolicy.AtomMapping getAtomMapping() {
        return atomMapping;
    }

    public static MatchingBehaviorStep of(GroupContainer reference, GroupContainer query) {
        return new MatchingBehaviorStep(reference, query);
    }

    public static class MatchingBehaviorStep {
        final GroupContainer reference;
        final GroupContainer query;
        AlignmentPolicy.AtomMapping atomMapping;

        MatchingBehaviorStep(GroupContainer reference, GroupContainer query) {
            this.reference = reference;
            this.query = query;
        }

        public StructureAlignmentQuery matchingBehavior(AlignmentPolicy.AtomMapping atomMapping) {
            this.atomMapping = atomMapping;
            return new StructureAlignmentQuery(this);
        }
    }
}

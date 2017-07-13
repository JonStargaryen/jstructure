package de.bioforscher.jstructure.align;

import de.bioforscher.jstructure.model.structure.container.GroupContainer;

/**
 * Builder for rather complex structure alignments.
 * Created by S on 12.07.2017.
 */
public class StructureAlignmentBuilder {
    public static MatchingBehaviorStep builder(GroupContainer reference, GroupContainer query) {
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

        public ManipulationBehaviorStep matchingBehavior(AlignmentPolicy.AtomMapping atomMapping) {
            this.atomMapping = atomMapping;
            return new ManipulationBehaviorStep(this);
        }
    }

    public static class ManipulationBehaviorStep {
        final GroupContainer reference;
        final GroupContainer query;
        final AlignmentPolicy.AtomMapping atomMapping;
        AlignmentPolicy.ManipulationBehavior manipulationBehavior;

        ManipulationBehaviorStep(MatchingBehaviorStep matchingBehaviorStep) {
            this.reference = matchingBehaviorStep.reference;
            this.query = matchingBehaviorStep.query;
            this.atomMapping = matchingBehaviorStep.atomMapping;
        }

        public StructureAlignmentStep manipulationBehavior(AlignmentPolicy.ManipulationBehavior manipulationBehavior) {
            this.manipulationBehavior = manipulationBehavior;
            return new StructureAlignmentStep(this);
        }
    }

    public static class StructureAlignmentStep {
        final GroupContainer reference;
        final GroupContainer query;
        final AlignmentPolicy.AtomMapping atomMapping;
        final AlignmentPolicy.ManipulationBehavior manipulationBehavior;

        StructureAlignmentStep(ManipulationBehaviorStep manipulationBehaviorStep) {
            this.reference = manipulationBehaviorStep.reference;
            this.query = manipulationBehaviorStep.query;
            this.atomMapping = manipulationBehaviorStep.atomMapping;
            this.manipulationBehavior = manipulationBehaviorStep.manipulationBehavior;
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

        public AlignmentPolicy.ManipulationBehavior getManipulationBehavior() {
            return manipulationBehavior;
        }
    }
}

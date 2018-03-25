package de.bioforscher.jstructure.align.query;

import de.bioforscher.jstructure.model.structure.Structure;

public class StructureAlignmentQuery implements AlignmentQuery {
    private final Structure reference;
    private final Structure query;

    private StructureAlignmentQuery(BuildStep buildStep) {
        this.reference = buildStep.reference;
        this.query = buildStep.query;
    }

    public Structure getReference() {
        return reference;
    }

    public Structure getQuery() {
        return query;
    }

    public static ReferenceStep builder() {
        return new ReferenceStep();
    }

    public static class ReferenceStep {
        public QueryStep reference(Structure reference) {
            return new QueryStep(reference);
        }
    }

    public static class QueryStep {
        private final Structure reference;

        QueryStep(Structure reference) {
            this.reference = reference;
        }

        public BuildStep query(Structure query) {
            return new BuildStep(reference, query);
        }
    }

    public static class BuildStep {
        private final Structure reference;
        private final Structure query;

        BuildStep(Structure reference, Structure query) {
            this.reference = reference;
            this.query = query;
        }

        public StructureAlignmentQuery build() {
            return new StructureAlignmentQuery(this);
        }
    }
}

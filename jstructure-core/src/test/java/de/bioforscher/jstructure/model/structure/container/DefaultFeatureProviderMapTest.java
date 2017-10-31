package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

public class DefaultFeatureProviderMapTest {
    private Structure structure;
    private Atom atom;

    @Before
    public void setup() {
        structure = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        atom = structure.atoms()
                .findFirst()
                .get();
    }

    @Test(expected = ComputationException.class)
    public void shouldFailForMissingAnnotation() {
        atom.getFeature(Feature1.class);
    }

    @Test
    public void shouldResolveFeatureProvider() {
        atom.getFeature(Feature2.class);
    }

    static class Feature1 extends FeatureContainerEntry {
        public Feature1() {
            super(null);
        }
    }

    @DefaultFeatureProvider(FeatureProvider2.class)
    static class Feature2 extends FeatureContainerEntry {
        public Feature2(FeatureProvider featureProvider) {
            super(featureProvider);
        }
    }

    static class FeatureProvider2 extends FeatureProvider {
        @Override
        protected void processInternally(Structure structure) {
            structure.atoms().forEach(atom -> atom.getFeatureContainer().addFeature(new Feature2(this)));
        }
    }
}
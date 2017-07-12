package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Checks especially the createCopy() methods of the data model.
 * Created by S on 23.11.2016.
 */
public class ModelIntegrityTest {
    private Protein protein;

    @Before
    public void setup() {
        this.protein = ProteinParser.source("1brr").parse();
    }

    @Test
    public void testGroupCreateCopy() {
        AtomContainer groupOriginal = protein.getGroups().get(0);
        AtomContainer groupCopy = groupOriginal.createCopy();
        Assert.assertTrue(groupCopy instanceof Group);
    }

    @Test
    public void testChainCreateCopy() {
        GroupContainer chainCopy = protein.getChains()
                .get(0)
                .createCopy();
        Assert.assertTrue(chainCopy instanceof Chain);
    }

    @Test
    public void testProteinCreateCopy() {
        ChainContainer proteinCopy = protein.createCopy();
        Assert.assertTrue(proteinCopy instanceof Protein);
    }

    @Test
    public void shouldNotCopyFeatureMapEntriesDuringCreateCopy() {
        protein.getFeatureContainer().addFeature(new AdditionalFeatureEntry(new AdditionalFeatureProvider()));
        Assert.assertTrue(protein.getFeatureContainer().getFeatureOptional(AdditionalFeatureEntry.class).isPresent());
        Assert.assertFalse(protein.createCopy().getFeatureContainer().getFeatureOptional(AdditionalFeatureEntry.class).isPresent());
    }

    public static class AdditionalFeatureEntry extends FeatureContainerEntry {
        public AdditionalFeatureEntry(AbstractFeatureProvider featureProvider) {
            super(featureProvider);
        }
    }

    @FeatureProvider(provides = AdditionalFeatureEntry.class)
    public static class AdditionalFeatureProvider extends AbstractFeatureProvider {
        @Override
        protected void processInternally(Protein protein) {

        }
    }
}

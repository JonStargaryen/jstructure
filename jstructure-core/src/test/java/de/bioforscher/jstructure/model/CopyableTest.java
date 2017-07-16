package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.Alanine;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Assert that the specification of createDeepCopy() and createShallowCopy() is correct.
 * Created by S on 16.07.2017.
 */
public class CopyableTest {
    private Structure protein;
    private Chain chain;
    private Group group;
    private Alanine alanine;
    private Atom atom;

    @Before
    public void setup() {
        protein = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        chain = protein.select()
                .asChain();
        group = chain.select()
                .asGroup();
        alanine = (Alanine) protein.select()
                .groupName("ALA")
                .asAminoAcid();
        atom = chain.select()
                .asAtom();
    }

    @Test
    public void shouldCloneConcreteImplementationsOfAminoAcids() {
        // clone containers led to groups losing their identity meaning they were no longer concrete amino acid
        // implementations but mere groups

        // will fail when copy does not contain amino acids
        ChainContainer proteinCopy = protein.createDeepCopy();
        System.out.println(proteinCopy.getAminoAcidSequence());
        Assert.assertFalse("no amino acid sequence - clone compromised", proteinCopy.getAminoAcidSequence().isEmpty());

        GroupContainer chainCopy = protein.getChains().get(0).createDeepCopy();
        System.out.println(chainCopy.getAminoAcidSequence());
        Assert.assertFalse("no amino acid sequence - clone compromised", chainCopy.getAminoAcidSequence().isEmpty());
    }

    @Test
    public void shouldCopyProteinShallow() {
        Structure copy = protein.createShallowCopy();
        // shallow copy should not contain references
        Assert.assertTrue(copy.getChains().isEmpty());
    }

    @Test
    public void shouldCopyProteinDeep() {
        Structure copy = protein.createDeepCopy();
        // deep copy should contain references, but those are also cloned
        Assert.assertFalse(copy.getChains().isEmpty());
        Chain clonedChain = copy.select().asChain();
        Assert.assertNotEquals(chain, clonedChain);
        Assert.assertNotEquals(chain.getParentStructure(), clonedChain.getParentStructure());
    }

    @Test
    public void shouldCopyChainShallow() {
        Chain copy = chain.createShallowCopy();
        // shallow copy should not contain references
        Assert.assertTrue(copy.getGroups().isEmpty());
    }

    @Test
    public void shouldCopyChainDeep() {
        Chain copy = chain.createDeepCopy();
        // deep copy should contain references, but those are also cloned
        Assert.assertFalse(copy.getGroups().isEmpty());
        Group clonedGroup = copy.select().asGroup();
        Assert.assertNotEquals(group, clonedGroup);
        Assert.assertNotEquals(group.getParentChain(), clonedGroup.getParentChain());
    }

    @Test
    public void shouldCopyAlanineShallow() {
        Alanine copy = (Alanine) alanine.createShallowCopy();
        // shallow copy should not contain references
        Assert.assertTrue(copy.getAtoms().isEmpty());
    }

    @Test
    public void shouldCopyAlanineDeep() {
        Alanine copy = (Alanine) alanine.createDeepCopy();
        // deep copy should contain references, but those are also cloned
        Assert.assertFalse(copy.getAtoms().isEmpty());
        Atom clonedAtom = copy.select().asAtom();
        Assert.assertNotEquals(atom, clonedAtom);
        Assert.assertNotEquals(atom.getParentGroup(), clonedAtom.getParentGroup());
    }

    @Test
    public void shouldCopyAtomShallow() {
        Atom copy = atom.createShallowCopy();
        // shallow copy should not contain references
        Assert.assertEquals(Group.UNKNOWN_GROUP, copy.getParentGroup());
    }

    @Test
    public void shouldCopyAtomDeep() {
        Atom copy = atom.createDeepCopy();
        // deep copy should contain references, but those to parent are not cloned
        Assert.assertEquals(atom.getParentGroup(), copy.getParentGroup());
    }

    @Test
    public void shouldNotCopyFeatureMapEntriesDuringCreateCopy() {
        protein.getFeatureContainer().addFeature(new AdditionalFeatureEntry(new AdditionalFeatureProvider()));
        Assert.assertTrue(protein.getFeatureContainer().getFeatureOptional(AdditionalFeatureEntry.class).isPresent());
        Assert.assertFalse(protein.createDeepCopy().getFeatureContainer().getFeatureOptional(AdditionalFeatureEntry.class).isPresent());
    }

    public static class AdditionalFeatureEntry extends FeatureContainerEntry {
        public AdditionalFeatureEntry(AbstractFeatureProvider featureProvider) {
            super(featureProvider);
        }
    }

    @FeatureProvider(provides = AdditionalFeatureEntry.class)
    public static class AdditionalFeatureProvider extends AbstractFeatureProvider {
        @Override
        protected void processInternally(Structure protein) {

        }
    }
}
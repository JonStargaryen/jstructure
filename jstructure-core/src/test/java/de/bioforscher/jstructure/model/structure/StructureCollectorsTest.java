package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Test capabilities of the structure collectors.
 * Created by bittrich on 5/29/17.
 */
public class StructureCollectorsTest {
    private Protein protein;

    @Before
    public void setup() {
        protein = ProteinParser.source("1brr").parse();
    }

    @Test
    public void shouldCloneConcreteImplementationsOfAminoAcids() {
        // clone containers led to groups losing their identity meaning they were no longer concrete amino acid
        // implementations but mere groups

        // will fail when copy does not contain amino acids
        ChainContainer proteinCopy = protein.createCopy();
        System.out.println(proteinCopy.getAminoAcidSequence());
        Assert.assertFalse("no amino acid sequence - clone compromised", proteinCopy.getAminoAcidSequence().isEmpty());

        GroupContainer chainCopy = protein.getChains().get(0).createCopy();
        System.out.println(chainCopy.getAminoAcidSequence());
        Assert.assertFalse("no amino acid sequence - clone compromised", chainCopy.getAminoAcidSequence().isEmpty());
    }

    @Test
    public void ensureOriginalInstanceIntegrityWhileClonesKeepParentReferences() {
        // cloning chains, groups or atoms should keep reference to the original instance without it being aware of cloned instances

        int initialAtomCount = protein.getAtoms().size();
        int initialGroupCount = protein.getGroups().size();

        List<Atom> clonedAtoms1 = protein.atoms()
                .map(Atom::new)
                .collect(Collectors.toList());

        AtomContainer collectedAtoms2 = protein.atoms()
                .collect(StructureCollectors.toAtomContainer());

        AtomContainer clonedAtoms3 = protein.getGroups().get(0).createCopy();

        // ensure number of instances associated to the original instances is intact
        Assert.assertEquals(initialAtomCount, protein.getAtoms().size());
        Assert.assertEquals(initialGroupCount, protein.getGroups().size());

        // ensure collected atoms point to the synthetic container

        // ensure that parent references up to the protein are intact though
        Assert.assertEquals(protein, clonedAtoms1.get(0).getParentGroup().getParentChain().getParentProtein());
        Assert.assertEquals(protein, collectedAtoms2.getAtoms().get(0).getParentGroup().getParentChain().getParentProtein());
        Assert.assertEquals(protein, clonedAtoms3.getAtoms().get(0).getParentGroup().getParentChain().getParentProtein());
    }

    @Test
    public void shouldNotManipulateStaticUnknownInstance() {
        // protein/chain/group contain 'unknown' instances which are the home to clone model instances without parent reference

        Group unknownGroup = Group.UNKNOWN_GROUP;
        Assert.assertTrue(unknownGroup.getAtoms().isEmpty());

        // collecting stuff to an AtomContainer will clone the unknown instance and assign the cloned atoms to that container
        AtomContainer atomContainer = protein.atoms()
                .collect(StructureCollectors.toAtomContainer());

        // references should be distinct
        Assert.assertNotEquals(unknownGroup, atomContainer);

        // original unknown group should still not contain any children
        Assert.assertTrue(unknownGroup.getAtoms().isEmpty());
    }

    @Test
    public void shouldCreateDistinctContainer() {
        // create copy
        AtomContainer original = protein.aminoAcids()
                .findFirst()
                .get();
        String initialPdbRecord = original.getPdbRepresentation();
        AtomContainer copy = original.createCopy();

        // manipulate coordinates of copy
        copy.getAtoms().get(0).setCoordinates(new double[] { 0, 0, 0 });

        // should not affect original
        Assert.assertNotEquals(initialPdbRecord, copy.getPdbRepresentation());
        Assert.assertEquals(initialPdbRecord, original.getPdbRepresentation());
    }

    @Test
    public void toAtomContainer() throws Exception {
        AtomContainer container = protein.atoms()
                .collect(StructureCollectors.toAtomContainer());
        container.atoms().forEach(atom -> Assert.assertEquals("parent reference was lost",
                "1brr",
                atom.getParentGroup().getParentChain().getParentProtein().getProteinIdentifier().getPdbId()));
    }

    @Test
    public void toGroupContainer() throws Exception {
        GroupContainer container = protein.groups()
                .collect(StructureCollectors.toGroupContainer());
        container.groups().forEach(group -> Assert.assertEquals("parent reference was lost",
                "1brr",
                group.getParentChain().getParentProtein().getProteinIdentifier().getPdbId()));
    }

    @Test
    public void toChainContainer() throws Exception {
        ChainContainer container = protein.chains()
                .collect(StructureCollectors.toChainContainer());
    }
}
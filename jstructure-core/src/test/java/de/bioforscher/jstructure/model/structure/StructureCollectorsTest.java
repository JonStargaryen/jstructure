package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

/**
 * Test capabilities of the structure collectors.
 * Created by bittrich on 5/29/17.
 */
public class StructureCollectorsTest {
    private Structure protein;

    @Before
    public void setup() {
        protein = StructureParser.fromInputStream(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1BRR))
                .minimalParsing(true)
                .parse();
    }

    @Test
    public void shouldNotManipulateStaticUnknownInstance() {
        // protein/chain/group contain 'unknown' instances which are the home to clone model instances without parent reference

        Group unknownGroup = Group.UNKNOWN_GROUP;
        Assert.assertTrue(unknownGroup.getAtoms().isEmpty());

        // collecting stuff to an AtomContainer will clone the unknown instance and assign the cloned atoms to that container
        AtomContainer atomContainer = protein.atoms()
                .collect(StructureCollectors.toIsolatedStructure());

        // references should be distinct
        Assert.assertNotEquals(unknownGroup, atomContainer);

        // original unknown group should still not contain any children
        Assert.assertTrue(unknownGroup.getAtoms().isEmpty());
    }

    @Test
    public void shouldCreateDistinctContainer() {
        // create copy
        AminoAcid original = protein.select()
                .aminoAcids()
                .asAminoAcid();
        String initialPdbRecord = original.getPdbRepresentation();
        AtomContainer copy = original.createDeepCopy();

        // manipulate coordinates of copy
        copy.getAtoms().get(0).setCoordinates(new double[] { 0, 0, 0 });

        // should not affect original
        Assert.assertNotEquals(initialPdbRecord, copy.getPdbRepresentation());
        Assert.assertEquals(initialPdbRecord, original.getPdbRepresentation());
    }

    @Test
    public void shouldKeepStructureReferenceAtom() throws Exception {
        AtomContainer container = protein.atoms()
                .collect(StructureCollectors.toIsolatedStructure());
        container.atoms().forEach(atom -> Assert.assertEquals("parent reference was lost",
                "1brr",
                atom.getParentGroup().getParentChain().getParentStructure().getProteinIdentifier().getPdbId()));
    }

    @Test
    public void shouldKeepStructureReferenceGroup() throws Exception {
        GroupContainer container = protein.groups()
                .collect(StructureCollectors.toIsolatedStructure());
        container.groups().forEach(group -> Assert.assertEquals("parent reference was lost",
                "1brr",
                group.getParentChain().getParentStructure().getProteinIdentifier().getPdbId()));
    }

    @Test
    public void shouldKeepStructureReferenceChain() throws Exception {
        Structure container = protein.chains()
                .collect(StructureCollectors.toIsolatedStructure());
        Assert.assertEquals("parent reference was lost",
                "1brr",
                container.getProteinIdentifier().getPdbId());
    }

    @Test
    public void shouldKeepOrdering() {
        Structure copy = protein.select()
                .aminoAcids()
                .asIsolatedStructure();
        List<Atom> atoms = copy.getAtoms();
        //TODO test flawed
        for(int i = 1; i < atoms.size(); i++) {
            Atom a1 = atoms.get(i - 1);
            Atom a2 = atoms.get(i);
            Assert.assertTrue("copy is unordered " + a2.getParentGroup() + " " +  a2 + " and " + a1.getParentGroup() + " " + a1,
                    a2.getPdbSerial() > a1.getPdbSerial());
        }
    }
}
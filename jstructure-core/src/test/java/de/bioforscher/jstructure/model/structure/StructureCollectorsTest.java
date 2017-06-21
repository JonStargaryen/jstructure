package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

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
    public void shouldCreateAminoAcids() {
        // clone containers led to groups losing their identity meaning they were no longer concrete amino acid
        // implementations but mere groups

        // will fail when copy does not contain amino acids
        ChainContainer copy = protein.createCopy();
        copy.getAminoAcidSequence();
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
        //TODO rework concept, has there to be some structure for atoms (i.e. organize selection again into chains e.g.)
        BasicAtomContainer container = protein.atoms()
                .collect(StructureCollectors.toAtomContainer());
        container.atoms().forEach(atom -> Assert.assertEquals("parent reference was lost",
                "1brr",
                atom.getParentGroup().getParentChain().getParentProtein().getPdbId().getPdbId()));
        Assert.assertEquals("parent reference does not match",
                container.getAtoms()
                .get(0)
                .getParentGroup()
                .getParentChain()
                .getParentProtein(),
                container.getOrigin());
    }

    @Test
    public void toGroupContainer() throws Exception {
        BasicGroupContainer container = protein.groups()
                .collect(StructureCollectors.toGroupContainer());
        container.groups().forEach(group -> Assert.assertEquals("parent reference was lost",
                "1brr",
                group.getParentChain().getParentProtein().getPdbId().getPdbId()));
        Assert.assertEquals("parent reference does not match",
                container.getGroups()
                        .get(0)
                        .getParentChain()
                        .getParentProtein(),
                container.getOrigin());
    }

    @Test
    public void toChainContainer() throws Exception {
        BasicChainContainer container = protein.chains()
                .collect(StructureCollectors.toChainContainer());
        container.chains().forEach(chain -> Assert.assertEquals("parent reference was lost",
                "1brr",
                chain.getParentProtein().getPdbId().getPdbId()));
        Assert.assertEquals("parent reference does not match",
                container.getChains()
                        .get(0)
                        .getParentProtein(),
                container.getOrigin());
    }
}
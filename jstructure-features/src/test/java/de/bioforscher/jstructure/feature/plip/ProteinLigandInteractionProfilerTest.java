package de.bioforscher.jstructure.feature.plip;

import de.bioforscher.jstructure.feature.plip.model.InteractionContainer;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class ProteinLigandInteractionProfilerTest {
    private ProteinLigandInteractionProfiler instance;

    @Before
    public void setup() {
        this.instance = ProteinLigandInteractionProfiler.getInstance();
    }

    @Test
    public void test_ligand_1eve() {
        Structure structure = StructureParser.fromPdbId("1eve").parse();
        Group ligand = structure.select()
                .chainId("A")
                .residueNumber(2001)
                .asGroup();
        Group res84 = structure.select().residueNumber(84).asGroup();
        Group res279 = structure.select().residueNumber(279).asGroup();
        Group res330 = structure.select().residueNumber(330).asGroup();

        InteractionContainer interactions = instance.annotateLigandInteractions(structure).getSubsetOfInteractions(ligand);
        Assert.assertFalse(interactions.getSubsetOfInteractions(res84).getPiStackingInteractions().isEmpty());
        Assert.assertFalse(interactions.getSubsetOfInteractions(res279).getPiStackingInteractions().isEmpty());
        Assert.assertFalse(interactions.getSubsetOfInteractions(res330).getPiCationInteractions().isEmpty());
    }

    @Test
    public void test_ligand_1h2t() {
        Structure structure = StructureParser.fromPdbId("1h2t").parse();
        Group ligand = structure.select()
                .chainId("Z")
                .residueNumber(1151)
                .asGroup();

        InteractionContainer interactions = instance.annotateLigandInteractions(structure).getSubsetOfInteractions(ligand);
        Assert.assertFalse(interactions.getInteractions().isEmpty());

        //TODO offset in renumbering of residues
//        Group res20 = structure.getGroups().get(21);
//        Group res43 = structure.select().residueNumber(43).asGroup();
//        Group res112 = structure.select().residueNumber(112).asGroup();
//        Group res116 = structure.select().residueNumber(116).asGroup();
//        Assert.assertFalse(interactions.getSubsetOfInteractions(res20).getPiStackingInteractions().isEmpty());
//        Assert.assertFalse(interactions.getSubsetOfInteractions(res43).getPiStackingInteractions().isEmpty());
//        Assert.assertFalse(interactions.getSubsetOfInteractions(res112).getHydrogenBonds().isEmpty());
//        Assert.assertFalse(interactions.getSubsetOfInteractions(res116).getSaltBridges().isEmpty());
    }
}
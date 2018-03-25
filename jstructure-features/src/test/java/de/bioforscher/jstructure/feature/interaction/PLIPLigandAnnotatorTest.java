package de.bioforscher.jstructure.feature.interaction;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Test;

public class PLIPLigandAnnotatorTest {
    @Test
    public void shouldAnnotateLigandInteractions() {
        Structure structure = StructureParser.fromInputStream(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ)).parse();
        PLIPLigandAnnotator plipLigandAnnotator = new PLIPLigandAnnotator();
        plipLigandAnnotator.process(structure);

        PLIPInteractionContainer plipInteractionContainer = structure.getChains()
                .get(0)
                .getFeature(PLIPInteractionContainer.class);

        Assert.assertFalse("did not annotate ligand interactions correctly",
                plipInteractionContainer.getPiStackings().isEmpty());
        Assert.assertFalse("did not annotate ligand interactions correctly",
                plipInteractionContainer.getWaterBridges().isEmpty());

        Assert.assertFalse("did not associate interactions correctly",
                structure.select()
                        .chainName("A")
                        .residueNumber(84)
                        .asAminoAcid()
                        .getFeature(PLIPInteractionContainer.class)
                        .getPiStackings()
                        .isEmpty());
    }
}
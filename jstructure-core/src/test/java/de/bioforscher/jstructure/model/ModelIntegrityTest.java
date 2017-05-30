package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.aminoacid.Asparagine;
import de.bioforscher.jstructure.model.structure.aminoacid.Histidine;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
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
    public void shouldGetGroupCopy() {
        GroupContainer copiedGroups = protein.select()
                .chainName("A", "C")
                .groupName(Asparagine.THREE_LETTER_CODE)
                .asGroupContainer()
                .createCopy();
        Assert.assertTrue(copiedGroups instanceof Chain);
        System.out.println(copiedGroups.getPdbRepresentation());
    }

    @Test
    public void shouldGetAtomCopy() {
        AtomContainer clonedSelectedAtoms = protein.select()
                .chainName("A", "B")
                .aminoAcids()
                .groupName(Histidine.THREE_LETTER_CODE)
                .atomName(AminoAcid.ALPHA_CARBON_NAME)
                .asAtomContainer()
                .createCopy();
        System.out.println(clonedSelectedAtoms.getPdbRepresentation());
    }

    @Test
    public void shouldGetProteinCopy() {
        ChainContainer copiedProtein = protein.createCopy();
        Assert.assertTrue(copiedProtein instanceof Protein);
    }
}

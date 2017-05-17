package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Checks especially the createCopy() methods of the data de.bioforscher.explorer.helices.model
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
        GroupContainer copiedGroups = Selection.on(protein)
                .chainName("A", "C")
                .aminoAcids(AminoAcidFamily.ASPARAGINE)
                .asGroupContainer()
                .createCopy();
        Assert.assertTrue(copiedGroups instanceof Chain);
        System.out.println(copiedGroups.composePDBRecord());
    }

    @Test
    public void shouldGetAtomCopy() {
        AtomContainer clonedSelectedAtoms = Selection.on(protein)
                .chainName("A", "B")
                .aminoAcids()
                .aminoAcids(AminoAcidFamily.HISTIDINE)
                .atomName(AminoAcidFamily.ATOM_NAMES.CA_ATOM_NAME)
                .asAtomContainer()
                .createCopy();
        System.out.println(clonedSelectedAtoms.composePDBRecord());
    }

    @Test
    public void shouldGetProteinCopy() {
        ChainContainer copiedProtein = protein.createCopy();
        Assert.assertTrue(copiedProtein instanceof Protein);
    }
}

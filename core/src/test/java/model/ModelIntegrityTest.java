package model;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Checks especially the getCopy() methods of the data model
 * Created by S on 23.11.2016.
 */
public class ModelIntegrityTest {
    private Protein protein;

    @Before
    public void setup() {
        this.protein = ProteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldGetGroupCopy() {
        GroupContainer copiedGroups = Selection.on(protein)
                .chainName("A", "C")
                .aminoAcids(AminoAcid.ASPARAGINE)
                .asGroupContainer()
                .getCopy();
        Assert.assertTrue(copiedGroups instanceof Group);
        System.out.println(copiedGroups.composePDBRecord());
    }

    @Test
    public void shouldGetAtomCopy() {
        AtomContainer clonedSelectedAtoms = Selection.on(protein)
                .chainName("A", "B")
                .aminoAcids()
                .aminoAcids(AminoAcid.HISTIDINE)
                .atomName(AminoAcid.ATOM_NAMES.CA_ATOM_NAME)
                .asAtomContainer()
                .getCopy();
        System.out.println(clonedSelectedAtoms.composePDBRecord());
    }

    @Test
    public void shouldGetProteinCopy() {
        ChainContainer copiedProtein = protein.getCopy();
        Assert.assertTrue(copiedProtein instanceof Protein);
    }
}

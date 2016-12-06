package selection;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Chain;
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
 * Validates the behaviour of the selection API.
 * Created by S on 21.11.2016.
 */
public class SelectionIntegrationTest {
    private Protein protein;

    @Before
    public void setup() {
        protein = ProteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldSelectChains() {
        System.out.println("Chain 'A':");
        Chain chainA = Selection.on(protein)
                .asChain();
        System.out.println(chainA);

        System.out.println("Chain 'B':");
        Chain chainB = Selection.on(protein)
                .chainName("B")
                .asChain();
        System.out.println(chainB);

        System.out.println("Chains 'B', 'C':");
        ChainContainer chains = Selection.on(protein)
                .chainName("B", "C")
                .asChainContainer();
        Assert.assertEquals(2, chains.getChains().size());

        chains.chains().forEach(System.out::println);
    }

    @Test
    public void shouldSelectGroupsFromChains() {
        // select all arginines in chain 'B' and 'C'
        System.out.println("arginines in Chain 'B' and 'C'");
        GroupContainer argGroups = Selection.on(protein)
                .chainName("B", "C")
                .aminoAcids(AminoAcid.ARGININE)
                .asGroupContainer();
        System.out.println(argGroups.composePDBRecord());

        // select hetatms in chain 'A'
        System.out.println("Hetatms in Chain 'B':");
        GroupContainer hets = Selection.on(protein)
                .chainName("C")
                .hetatms()
                .asGroupContainer();
        // assert that any atoms where selected
        Assert.assertTrue(hets.getAtoms().size() > 0);
        System.out.println(hets.composePDBRecord());
    }

    @Test
    public void shouldSelectAtomsFromGroupsFromChains() {
        System.out.println("alpha carbons and beta carbons of PHEs in chain 'A' and 'C'");
        AtomContainer atoms = Selection.on(protein)
                .chainName("A", "C")
                .aminoAcids(AminoAcid.PHENYLALANINE)
                .atomName(AminoAcid.ATOM_NAMES.CA_ATOM_NAME, AminoAcid.ATOM_NAMES.CB_ATOM_NAME)
                .asAtomContainer();
        System.out.println(atoms.composePDBRecord());
    }
}

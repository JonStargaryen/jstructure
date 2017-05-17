package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Validates the behaviour of the selection API.
 * Created by S on 21.11.2016.
 */
public class SelectionTest {
    private Protein protein;

    @Before
    public void setup() {
        protein = ProteinParser.source("1brr").parse();
    }

    @Test
    public void shouldSelectBasedOnRange() {
        List<Group> selection = Selection.on(protein)
                .chainName("A")
                .residueNumber(new IntegerRange(100, 105))
                .asFilteredGroups()
                .collect(Collectors.toList());

        selection.forEach(group -> Assert.assertTrue(group.getResidueNumber() >= 100 && group.getResidueNumber() <= 105));
        System.out.println(selection);
        Assert.assertTrue(selection.size() == 6);
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
                .aminoAcids(AminoAcidFamily.ARGININE)
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
                .aminoAcids(AminoAcidFamily.PHENYLALANINE)
                .atomName(AminoAcidFamily.ATOM_NAMES.CA_ATOM_NAME, AminoAcidFamily.ATOM_NAMES.CB_ATOM_NAME)
                .asAtomContainer();
        System.out.println(atoms.composePDBRecord());
    }

    private static final String givenName = "ThereIsAlwaysMoneyInTheBananaStand";

    @Test
    public void shouldNameContainer() {
        GroupContainer container = Selection.on(protein)
                .aminoAcids()
                .nameContainer(givenName)
                .asGroupContainer();

        Assert.assertTrue(container.getIdentifier().equals(givenName));
    }

    @Test
    public void shouldNameAfterParentContainer() {
        GroupContainer parentContainer = Selection.on(protein)
                .cloneElements()
                .nameContainer(givenName)
                .asGroupContainer();

        GroupContainer container = Selection.on(parentContainer)
                .aminoAcids()
                .asGroupContainer();

        Assert.assertTrue(container.getIdentifier().equals(givenName));
    }

    @Test
    public void shouldFindSimilarNeighbors() {
        Protein protein = ProteinParser.source("1acj").parse();

        Atom firstAtom = protein.getAtoms().get(0);
        double probeDistance = 4.4499;

        int neighboringGroupCountSelectionAPI = (int) Selection.on(protein)
                .atomSelection()
                .distance(firstAtom, probeDistance)
                .asFilteredAtoms()
                .peek(System.out::println)
                .count();

        int neighboringGroupCountNaive = (int) protein.atoms()
                .filter(atom -> LinearAlgebra3D.distance(atom.getCoordinates(), firstAtom.getCoordinates()) < probeDistance)
                .peek(System.out::println)
                .count();

        Assert.assertEquals(neighboringGroupCountNaive, neighboringGroupCountSelectionAPI);
    }
}

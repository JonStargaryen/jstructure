package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.aminoacid.Arginine;
import de.bioforscher.jstructure.model.structure.aminoacid.Phenylalanine;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.Optional;
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
    public void shouldSelectAminoAcidsAndReturnAsContainer() {
        List<Group> groups = protein.select()
                .aminoAcids()
                .asFilteredGroups()
                .collect(Collectors.toList());

        groups.forEach(group -> Assert.assertTrue("group " + group + " of original selected stream is no amino acid, was " + group.getClass().getSimpleName(), group instanceof AminoAcid));

        GroupContainer container = groups.stream()
                .collect(StructureCollectors.toGroupContainer());
        container.groups().forEach(group -> Assert.assertTrue("group " + group + " of container stream is no amino acid, was " + group.getClass().getSimpleName(), group instanceof AminoAcid));
    }

    @Test
    public void shouldRetainParentReferences() {
        Optional<Group> group = protein.select()
                .chainName("A")
                .residueNumber(100)
                .asOptionalGroup();

        System.out.println(group.get().getParentChain().getParentProtein().getIdentifier());
    }

    @Test
    public void shouldSelectBasedOnRange() {
        List<Group> selection = protein.select()
                .chainName("A")
                .residueNumber(new IntegerRange(100, 105))
                .asFilteredGroups()
                .collect(Collectors.toList());

        selection.forEach(group -> Assert.assertTrue(group.getResidueNumber().getResidueNumber() >= 100 &&
                group.getResidueNumber().getResidueNumber() <= 105));
        System.out.println(selection);
        Assert.assertTrue(selection.size() == 6);
    }

    @Test
    public void shouldSelectChains() {
        System.out.println("Chain 'A':");
        Chain chainA = protein.select()
                .asChain();
        System.out.println(chainA);

        System.out.println("Chain 'B':");
        Chain chainB = protein.select()
                .chainName("B")
                .asChain();
        System.out.println(chainB);

        System.out.println("Chains 'B', 'C':");
        List<Chain> chains = protein.select()
                .chainName("B", "C")
                .asFilteredChains()
                .collect(Collectors.toList());
        Assert.assertEquals("chains do not match, expected 'B' and 'C' - found: "  + chains, 2, chains.size());

        chains.forEach(System.out::println);
    }

    @Test
    public void shouldSelectGroupsFromChains() {
        // select all arginines in chain 'B' and 'C'
        System.out.println("arginines in Chain 'B' and 'C'");
        GroupContainer argGroups = protein.select()
                .chainName("B", "C")
                .groupName(Arginine.THREE_LETTER_CODE)
                .asGroupContainer();
        System.out.println(argGroups.getPdbRepresentation());

        // select hetatms in chain 'A'
        System.out.println("Hetatms in Chain 'B':");
        GroupContainer hets = protein.select()
                .chainName("C")
                .hetatms()
                .asGroupContainer();
        // assert that any atoms where selected
        Assert.assertTrue(hets.getAtoms().size() > 0);
        System.out.println(hets.getPdbRepresentation());
    }

    @Test
    public void shouldSelectAtomsFromGroupsFromChains() {
        System.out.println("alpha carbons and beta carbons of PHEs in chain 'A' and 'C'");
        AtomContainer atoms = protein.select()
                .chainName("A", "C")
                .groupName(Phenylalanine.THREE_LETTER_CODE)
                .atomName(AminoAcid.ALPHA_CARBON_NAME, AminoAcid.BETA_CARBON_NAME)
                .asAtomContainer();
        System.out.println(atoms.getPdbRepresentation());
    }

    private static final String givenName = "ThereIsAlwaysMoneyInTheBananaStand";

    @Test
    public void shouldNameContainer() {
        GroupContainer container = protein.select()
                .aminoAcids()
                .nameContainer(givenName)
                .asGroupContainer();

        Assert.assertTrue(container.getIdentifier().equals(givenName));
    }

    @Test
    public void shouldNameAfterParentContainer() {
        GroupContainer parentContainer = protein.select()
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

        Atom atom = protein.getAtoms().get(100);
        double probeDistance = 4;

        List<Atom> neighboringGroupCountSelectionAPI = protein.select()
                .atomSelection()
                .atomDistance(atom, probeDistance)
                .asFilteredAtoms()
                .collect(Collectors.toList());

        List<Atom> neighboringGroupCountNaive = protein.atoms()
                .filter(a -> a.calculate().distance(atom.getCoordinates()) < probeDistance)
                .collect(Collectors.toList());

        Assert.assertTrue(neighboringGroupCountSelectionAPI.containsAll(neighboringGroupCountNaive));
        Assert.assertTrue(neighboringGroupCountNaive.containsAll(neighboringGroupCountSelectionAPI));
    }
}

package reconstruction;

import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.reconstruction.sidechain.pulchra.PULCHRA;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Checks PULCHRA implementation.
 * Created by S on 04.11.2016.
 */
public class PULCHRAFunctionalTest {
    private Protein protein;
    private Protein proteinCopy;

    @Before
    public void setup() {
        // for sake of simplicity a structure with a single getChain
        protein = ProteinParser.source("1acj").parse();
        proteinCopy = ProteinParser.source("1acj").parse();
    }

    @Test
    public void shouldMutateGivenAminoAcidPositionToGlycine() {
        Protein protein = ProteinParser.source("1brr").parse();
        Group group = protein.select()
                .chainName("A")
                .residueNumber(100)
                .asGroup();
        System.out.println("original entry:" + System.lineSeparator() + group.composePDBRecord());
        new PULCHRA().mutatePosition(group, "GLY");
        System.out.println("mutated entry:" + System.lineSeparator() + group.composePDBRecord());
    }

    @Test
    public void shouldMutateGivenAminoAcidPositionToTryptophane() {
        Protein protein = ProteinParser.source("1brr").parse();
        Group group = protein.select()
                .chainName("A")
                .residueNumber(100)
                .asGroup();
        System.out.println("original entry:" + System.lineSeparator() + group.composePDBRecord());
        new PULCHRA().mutatePosition(group, "TRP");
        System.out.println("mutated entry:" + System.lineSeparator() + group.composePDBRecord());
    }

    @Test
    public void shouldCreatePULCHRAInstance() {
        new PULCHRA();
    }

    @Test
    public void checkIfRebuildIsComplete() {
        // initial count
        int proteinAtomCount = determineAtomCount(protein);

        // forget everything but alpha carbons
        GroupContainer residues = Selection.on(proteinCopy)
                .aminoAcids()
                .asGroupContainer();
        List<Atom> alphaCarbons = Selection.on(residues)
                .alphaCarbonAtoms()
                .asFilteredAtoms()
                .collect(Collectors.toList());
        proteinCopy.clearAtoms();
        Combinatorics.sequentialPairsOf(residues.getGroups(), alphaCarbons)
            .forEach(pair -> pair.getLeft().addAtom(pair.getRight()));

        // reconstruct and count again - all atoms should be reconstructed
        new PULCHRA().reconstruct(proteinCopy);
        int reconstructedAtomCount = determineAtomCount(proteinCopy);

        // check for the atom name sequence to match
        Assert.assertEquals(determineAtomSequence(protein), determineAtomSequence(proteinCopy));

        // check for the number of atoms to match
        System.out.printf("count before reconsruction: %d - after: %d", proteinAtomCount, reconstructedAtomCount);
        Assert.assertEquals(proteinAtomCount, reconstructedAtomCount);
    }

    private String determineAtomSequence(Protein protein) {
        return protein.atoms()
                      .map(Atom::getName)
                      .collect(Collectors.joining());
    }

    private int determineAtomCount(Protein protein) {
        return protein.getAtoms().size();
    }
}

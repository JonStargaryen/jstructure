package reconstruction;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.reconstruction.sidechain.pulchra.PULCHRA;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.Optional;
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
        protein = ProteinParser.parseProteinById("1acj");
        proteinCopy = ProteinParser.parseProteinById("1acj");
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
        List<Residue> residues = proteinCopy.residues().collect(Collectors.toList());
        List<Atom> alphaCarbons =  residues.stream()
               .map(Residue::findAlphaCarbon)
               .map(Optional::get)
               .collect(Collectors.toList());
        proteinCopy.clear();
        Pair.sequentialPairsOf(residues, alphaCarbons)
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
        return (int) protein.residues()
                            .flatMap(Residue::atoms)
                            .count();
    }
}

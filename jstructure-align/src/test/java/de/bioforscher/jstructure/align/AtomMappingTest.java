package de.bioforscher.jstructure.align;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

/**
 * Tests whether atom mapping correctly pairs atoms.
 * Created by S on 16.07.2017.
 */
public class AtomMappingTest {
    private Structure protein1;
    private Structure protein2;

    @Before
    public void setup() {
        protein1 = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        protein2 = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
    }

    @Test
    public void testComparableAtomNames() {
        AlignmentPolicy.AtomMapping atomMapping = AlignmentPolicy.MatchingBehavior.comparableAtomNames;
        List<Pair<Atom, Atom>> mappedAtoms = atomMapping.determineAtomMapping(protein1, protein2);

        Assert.assertEquals("pairs where missing for identical proteins",
                protein1.atoms().count(),
                mappedAtoms.size());
        mappedAtoms
            .forEach(pair -> Assert.assertTrue("pairs were not of equal atom names",
                    pair.getLeft().getName().equals(pair.getRight().getName())));
    }
}
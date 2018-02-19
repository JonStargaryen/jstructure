package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.mathematics.Fragment;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests the combinatoric functions of Pair and Fragment.
 * Asserts are realized by mere size checks, at the moment no underlying integrity is checked.
 * Created by S on 02.10.2016.
 */
public class SetOperationsTest {
    private List<Integer> elements1;
    private List<Integer> elements2;

    @Before
    public void setup() {
        elements1 = IntStream.range(0, 11).boxed().collect(Collectors.toList());
        elements2 = IntStream.range(20, 31).boxed().collect(Collectors.toList());
    }

    @Test
    public void shouldFragmentCollection() {
        final int fragmentSize = 3;
        List<Fragment<Integer>> fragments = SetOperations.fragmentsOf(elements1, fragmentSize).collect(Collectors.toList());
        // fragmenting n elements into groups of fragmentSize should result in n - fragmentSize + 1 fragments
        Assert.assertEquals(fragments.size(), elements1.size() - fragmentSize + 1);
        fragments.forEach(System.out::println);
    }

    @Test
    public void shouldGenerateUnorderedPairs() {
        List<Pair<Integer, Integer>> pairs = SetOperations.unorderedPairsOf(elements1).collect(Collectors.toList());
        // generation of unordered pairs for n elements should return n * (n - 1) / 2
        Assert.assertEquals(pairs.size(), elements1.size() * (elements1.size() - 1) / 2);
        pairs.forEach(System.out::println);
    }

    @Test
    public void shouldComposeCartesianProduct() {
        List<Pair<Integer, Integer>> pairs = SetOperations.cartesianProductOf(elements1, elements2).collect(Collectors.toList());
        // cartesian product of N x M should return |N| * |M| elements
        Assert.assertEquals(pairs.size(), elements1.size() * elements2.size());
        pairs.forEach(System.out::println);
    }

    @Test
    public void testPairContains() {
        Double left1 = 10.0;
        String left2 = "true";
        SetOperationsTest right1 = this;
        Object right2 = new Object();

        Pair<Double, SetOperationsTest> pair1 = new Pair<>(left1, right1);
        Pair<String, Object> pair2 = new Pair<>(left2, right2);

        // should contain object no matter of order
        Assert.assertTrue(pair1.contains(left1));
        Assert.assertTrue(pair1.contains(right1));
        Assert.assertTrue(pair2.contains(left2));
        Assert.assertTrue(pair2.contains(right2));

        // should not contain other objects
        Assert.assertFalse(pair1.contains(left2));
        Assert.assertFalse(pair1.contains(right2));
        Assert.assertFalse(pair2.contains(left1));
        Assert.assertFalse(pair2.contains(right1));
    }
}
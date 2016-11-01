package mathematics;

import de.bioforscher.jstructure.model.Fragment;
import de.bioforscher.jstructure.model.Pair;
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
public class CombinatoricsFunctionalTest {
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
        List<Fragment<Integer>> fragments = Fragment.fragmentsOf(elements1, fragmentSize).collect(Collectors.toList());
        // fragmenting n elements into groups of fragmentSize should result in n - fragmentSize + 1 fragments
        Assert.assertEquals(fragments.size(), elements1.size() - fragmentSize + 1);
        fragments.forEach(System.out::println);
    }

    @Test
    public void shouldGenerateUnorderedPairs() {
        List<Pair<Integer, Integer>> pairs = Pair.uniquePairsOf(elements1).collect(Collectors.toList());
        // generation of unordered pairs for n elements should return n * (n - 1) / 2
        Assert.assertEquals(pairs.size(), elements1.size() * (elements1.size() - 1) / 2);
        pairs.forEach(System.out::println);
    }

    @Test
    public void shouldComposeCartesianProduct() {
        List<Pair<Integer, Integer>> pairs = Pair.cartesianProductOf(elements1, elements2).collect(Collectors.toList());
        // cartesian product of N x M should return |N| * |M| elements
        Assert.assertEquals(pairs.size(), elements1.size() * elements2.size());
        pairs.forEach(System.out::println);
    }
}

package combinatoric;

import de.bioforscher.jstructure.model.Fragment;
import de.bioforscher.jstructure.model.Pair;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO: this should be actual automatic tests
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
        Fragment.fragmentsOf(elements1, 3).forEach(System.out::println);
    }

    @Test
    public void shouldGenerateUnorderedPairs() {
        Pair.unorderedPairsOf(elements1).forEach(System.out::println);
    }

    @Test
    public void shouldComposeCartesianProduct() {
        Pair.cartesianProductOf(elements1, elements2).forEach(System.out::println);
    }
}

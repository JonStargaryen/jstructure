package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class ClusteringScorerTest {
    private List<Module> clustering1;
    private List<Module> clustering2;

    @Before
    public void setup() {
        this.clustering1 = createClusterFromModulesDat(TestUtils.getResourceAsLines("modularity/3cyt_I_generous.modules.dat"));
        // for definitions see Bai, 1995
        this.clustering2 = createClusterManually(new int[] { 70, 85 },
                new int[] { 36, 61 },
                new int[] { 20, 35 },
                new int[] { 62, 69 },
                new int[] { 1, 19 },
                new int[] { 86, 103 });
    }

    private List<Module> createClusterFromModulesDat(List<String> lines) {
        return lines.stream()
                .filter(line -> !line.startsWith("#"))
                .map(line -> line.split("---")[1])
                .map(String::trim)
                .map(line -> Pattern.compile("\\s+").splitAsStream(line)
                        .collect(Collectors.toList()))
                .map(Module::new)
                .collect(Collectors.toList());
    }

    private List<Module> createClusterManually(int[]... ranges) {
        return Stream.of(ranges)
                .map(range -> IntStream.range(range[0], range[1] + 1)
                        .mapToObj(String::valueOf)
                        .collect(Collectors.toList()))
                .map(Module::new)
                .collect(Collectors.toList());
    }

    @Test
    public void shouldComputeChiSquaredCoefficient() {
        System.out.println("chi-squared: " + ClusteringScorer.chiSquaredCoefficient(clustering1, clustering2));
    }

    @Test
    public void shouldComputeNaiveScore() {
        System.out.println("naive score: " + ClusteringScorer.naiveScore(clustering1, clustering2));
    }
}
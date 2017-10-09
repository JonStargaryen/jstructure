package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ClusteringScorerTest {
    private Pair<List<Module>, List<Module>> _1bf2;
    private Pair<List<Module>, List<Module>> _1f21;
    private Pair<List<Module>, List<Module>> _3cyt;
    private Pair<List<Module>, List<Module>> _3gkh;

    @Before
    public void setup() {
        this._1bf2 = new Pair<>(createClusterManually("A:179-180,476-480,506-506,559-572,576-576,578-578,582-582,615-748",
                        "B:181-183,327-365,372-372,374-432,434-455,489-502",
                        "C:184-326,367-371,433-433,580-580,602-603",
                        "D:456-475,481-488,503-505,507-516,573-575,577-577,579-579,581-581,604-606"),
                createClusterFromModulesDat(TestUtils.getResourceAsLines("modularity/1bf2_A_plip.modules.dat")));
        this._1f21 = new Pair<>(createClusterManually("A:5-13,18-28,31-42,127-143",
                        "B:43-58,64-69",
                        "C:71-88",
                        "D:100-112,115-120"),
                createClusterFromModulesDat(TestUtils.getResourceAsLines("modularity/1f21_A_plip.modules.dat")));
        this._3cyt = new Pair<>(createClusterManually("A:70-85",
                        "B:36-61",
                        "C:20-35,62-69",
                        "D:1-19,86-103"),
                createClusterFromModulesDat(TestUtils.getResourceAsLines("modularity/3cyt_I_plip.modules.dat")));
        this._3gkh = new Pair<>(createClusterManually("A:23-31,41-42,44-85,109-109,112-119,138-138",
                        "B:32-40,43-43,86-92,125-137,139-142,196-196",
                        "C:197-202",
                        "D:189-189,193-193,224-246",
                        "E:93-106,155-188,191-192,195-195",
                        "F:107-108,110-111,120-124,143-154,190-190,194-194,203-223"),
                createClusterFromModulesDat(TestUtils.getResourceAsLines("modularity/3gkh_A_plip.modules.dat")));
    }

    private List<Module> createClusterFromModulesDat(List<String> lines) {
        return lines.stream()
                .filter(line -> !line.startsWith("#"))
                .map(line -> new Module(line.split("\\s+")[0],
                        Pattern.compile("\\s+").splitAsStream(line.split("---")[1].trim())
                                .collect(Collectors.toList())))
                .collect(Collectors.toList());
    }

    private List<Module> createClusterManually(String... ranges) {
        List<Module> modules = new ArrayList<>();
        for(String range : ranges) {
            String id = range.split(":")[0];
            String rawRanges = range.split(":")[1];

            List<String> nodes = Pattern.compile(",").splitAsStream(rawRanges)
                    .map(rawRange -> rawRange.split("-"))
                    .flatMap(rawRange -> IntStream.range(Integer.valueOf(rawRange[0]), Integer.valueOf(rawRange[1]) + 1)
                            .mapToObj(String::valueOf))
                    .collect(Collectors.toList());

            modules.add(new Module(id, nodes));
        }

        return modules;
    }

    @Test
    public void shouldComputeChiSquaredCoefficient() {
        System.out.println("[1bf2] chi-squared: " + ClusteringScorer.chiSquaredCoefficient(_1bf2.getLeft(), _1bf2.getRight()));
        System.out.println("[1f21] chi-squared: " + ClusteringScorer.chiSquaredCoefficient(_1f21.getLeft(), _1f21.getRight()));
        System.out.println("[3cyt] chi-squared: " + ClusteringScorer.chiSquaredCoefficient(_3cyt.getLeft(), _3cyt.getRight()));
        System.out.println("[3gkh] chi-squared: " + ClusteringScorer.chiSquaredCoefficient(_3gkh.getLeft(), _3gkh.getRight()));
    }

    @Test
    public void shouldComputeNaiveScore() {
        System.out.println("[1bf2] naive score: " + ClusteringScorer.naiveScore(_1bf2.getLeft(), _1bf2.getRight()));
        System.out.println("[1f21] naive score: " + ClusteringScorer.naiveScore(_1f21.getLeft(), _1f21.getRight()));
        System.out.println("[3cyt] naive score: " + ClusteringScorer.naiveScore(_3cyt.getLeft(), _3cyt.getRight()));
        System.out.println("[3gkh] naive score: " + ClusteringScorer.naiveScore(_3gkh.getLeft(), _3gkh.getRight()));
    }

    @Test
    public void shouldComputeRandIndex() {
        System.out.println("[1bf2] rand-index: " + ClusteringScorer.randIndex(_1bf2.getLeft(), _1bf2.getRight()));
        System.out.println("[1f21] rand-index: " + ClusteringScorer.randIndex(_1f21.getLeft(), _1f21.getRight()));
        System.out.println("[3cyt] rand-index: " + ClusteringScorer.randIndex(_3cyt.getLeft(), _3cyt.getRight()));
        System.out.println("[3gkh] rand-index: " + ClusteringScorer.randIndex(_3gkh.getLeft(), _3gkh.getRight()));
    }
}
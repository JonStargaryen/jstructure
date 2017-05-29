package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The functional test for the energy profile calculator.
 * Created by bittrich on 12/15/16.
 */
public class EnergyProfileCalculatorTest {
    private List<String> ids;
    private DecimalFormat decimalFormat;
    private AbstractFeatureProvider featureProvider;

    @Before
    public void setup() {
        ids = Stream.of("1bs2", "1f7u", "1f7v", "3fnr", "4oby", "4q2t", "4q2x", "4q2y", "4r3z", "1li5", "1li7", "1u0b", "3c8z", "3sp1", "3tqo", "1euq", "1euy", "1exd", "1gsg", "1gtr", "1gts", "1nyl", "1o0b", "1qrs", "1qrt", "1qru", "1qtq", "1zjw", "2hz7", "2rd2", "2re8", "4h3s", "4jxx", "4jxz", "4jyz", "4p2b", "4r3z", "4ye6", "4ye8", "4ye9", "5bnz", "1g59", "1gln", "1j09", "1n75", "1n77", "1n78", "1nzj", "2cfo", "2cuz", "2cv0", "2cv1", "2cv2", "2dxi", "2ja2", "2o5r", "3afh", "3aii", "3akz", "3pnv", "3pny", "4gri", "1ile", "1jzq", "1jzs", "3zgz", "3ziu", "3zjt", "3zju", "3zjv", "4aq7", "4arc", "4ari", "4as1", "4cqn", "5ah5", "1a8h", "1f4l", "1p7p", "1pfv", "1pfw", "1pfy", "1pg0", "1pg2", "1qqt", "1woy", "2csx", "2ct8", "2d54", "2d5b", "2x1l", "2x1m", "3h97", "3h99", "3h9b", "3kfl", "3vu8", "4dlp", "4eg1", "4eg3", "4eg4", "4eg5", "4eg6", "4eg7", "4eg8", "4ega", "4mvw", "4mvx", "4mvy", "4mw0", "4mw1", "4mw2", "4mw4", "4mw5", "4mw6", "4mw7", "4mw9", "4mwb", "4mwc", "4mwd", "4mwe", "4py2", "1o5t", "1r6t", "1r6u", "1ulh", "1yi8", "1yia", "1yid", "2a4m", "2ake", "2azx", "2dr2", "2el7", "2g36", "2ip1", "2quh", "2qui", "2quj", "2quk", "2yy5", "3a04", "3a05", "3foc", "3hv0", "3i05", "3kt0", "3kt3", "3kt6", "3kt8", "3m5w", "3n9i", "3prh", "3sz3", "3tze", "3tzl", "4j75", "4j76", "4jfa", "1j1u", "1n3l", "1q11", "1u7d", "1u7x", "1vbm", "1vbn", "1wq3", "1wq4", "1x8x", "2cya", "2cyb", "2dlc", "2j5b", "2jan", "2pid", "2rkj", "2yxn", "2zp1", "3p0h", "3p0i", "3p0j", "3vgj", "3zxi", "4hjr", "4hjx", "4hk4", "4hpw", "4nd6", "4nd7", "4nda", "4ojm", "4q93", "4qbt", "1riq", "1yfr", "1yfs", "1yft", "1ygb", "2ztg", "3htz", "3hxu", "3hxv", "3hxw", "3hxx", "3hxy", "3hxz", "3hy0", "3hy1", "3wqy", "3wqz", "4xem", "4xeo", "1nnh", "3m4p", "3m4q", "1asy", "1asz", "1b8a", "1c0a", "1eov", "1eqr", "1il2", "1n9w", "1wyd", "3i7f", "3nel", "3nem", "3nen", "4ah6", "4j15", "4o2d", "4rmf", "1ati", "1b76", "1ggm", "2pme", "2pmf", "2q5h", "2q5i", "2zt5", "2zt6", "2zt7", "2zt8", "2zxf", "3rgl", "3ufg", "4h2x", "4kqe", "4kr2", "4kr3", "4qei", "1htt", "1kmm", "1kmn", "2el9", "1bbu", "1bbw", "1e1o", "1e1t", "1e22", "1e24", "1lyl", "3a5y", "3a5z", "3bju", "3g1z", "4dpg", "4h02", "4pg3", "4up7", "4up8", "4up9", "4upa", "4ycu", "4ycv", "4ycw", "1eiy", "1jjc", "1pys", "3cmq", "3hfv", "3hfz", "3l4g", "3pco", "3teg", "3tup", "4p71", "4p72", "4p73", "4p74", "4p75", "4tva", "1h4q", "1h4s", "1h4t", "1hc7", "1nj1", "1nj2", "1nj5", "1nj6", "1nj8", "2i4l", "2i4m", "2i4n", "2i4o", "2j3l", "2j3m", "3ial", "4hvc", "4k86", "4k87", "4k88", "4ncx", "4olf", "4twa", "4ydq", "2e3c", "2q7e", "2q7g", "2q7h", "2zce", "2zim", "2zin", "2zio", "3qtc", "3vqv", "3vqw", "3vqx", "4bw9", "4bwa", "4ch3", "4ch4", "4ch5", "4ch6", "4cs2", "4cs3", "4cs4", "4q6g", "4tqd", "4tqf", "2du3", "2du4", "2du5", "2du6", "2du7", "2odr", "1wle", "2dq3", "3lsq", "3lss", "3mey", "3mf1", "3mf2", "3qne", "3qo5", "3qo7", "3qo8", "3vbb", "3w3s", "4l87", "4rqe", "4rqf", "1evk", "1evl", "1fyf", "1kog", "1qf6", "3a31", "3a32", "3ugq", "3ugt", "3uh0", "4eo4", "4hwo", "4hwp", "4hwr", "4hws", "4hwt", "4p3n", "4p3o", "4p3p", "4ttv", "4yye")
               .collect(Collectors.toList());
        decimalFormat = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.US));
        featureProvider = FeatureProviderRegistry.resolve(EnergyProfile.class);
    }

    @Test
    @Ignore("skipping long calculation")
    public void shouldComputeEnergyProfilesForAllStructures() throws IOException {
        System.err.println("skipping energy profile calculation on aaRS data set");
        // this takes ages - uncomment if needed
        ids.stream()
                .map(path -> ProteinParser.source(path).parse())
                .peek(System.out::println)
                .forEach(protein -> {
                    featureProvider.process(protein);

                    System.out.println(protein.aminoAcids()
                            .map(Group::getFeatureContainer)
                            .map(container -> container.getFeature(EnergyProfile.class))
                            .mapToDouble(EnergyProfile::getSolvationEnergy)
                            .mapToObj(decimalFormat::format)
                            .collect(Collectors.joining(", ", "energy values: ", "")));
                });
    }

    @Test
    public void shouldProcessStructureWithSelenomethionine() {
        Protein protein = ProteinParser.source("3TQO").parse();
        featureProvider.process(protein);
    }

    @Test
    public void shouldProcessStructureWithMPH() {
        Protein protein = ProteinParser.source("1P7P").parse();
        featureProvider.process(protein);
    }

    @Test
    public void shouldProcessStructure() {
        Protein protein = ProteinParser.source("1ATI").parse();
        featureProvider.process(protein);
    }

    @Test
    public void shouldAgreeRegardingEnergyTerms() throws IOException {
        Protein protein = ProteinParser.source("1bs2").parse();
        featureProvider.process(protein);

        Files.lines(Paths.get(TestUtils.getResourceAsFilepath("energy/1bs2.ep2")))
                .filter(line -> line.startsWith("ENGY"))
                .map(line -> line.split("\t"))
                .forEach(split -> {
                    Group group = protein.select()
                            .chainName(split[1])
                            .residueNumber(Integer.valueOf(split[2]))
                            .asGroup();

                    Assert.assertEquals("energy values differ for " + group,
                            Double.valueOf(split[5]),
                            group.getFeatureContainer().getFeature(EnergyProfile.class).getSolvationEnergy(),
                            0.001);
                });
    }
}
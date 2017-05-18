package studies.gmlvq;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import studies.StudyConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Create some data set for the GMLVQ-project.
 * Created by bittrich on 4/11/17.
 */
public class ComposeItemSetMinerDataSet {
    private static final boolean FUNCTIONAL = false;
    private static final String FILENAME = "studies/gmlvq/" + (FUNCTIONAL ? "positives.txt" : "negatives.txt");

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####");
    private static final List<AbstractFeatureProvider> featureProviders = Stream.of(AccessibleSurfaceArea.class,
            LoopFraction.class,
            EnergyProfile.class)
            .map(FeatureProviderRegistry::resolve)
            .collect(Collectors.toList());

    public static void main(String[] args) throws IOException {
        int numberOfGroups = Files.lines(Paths.get(StudyConstants.getResourceAsFilepath(FILENAME))).findFirst().orElseThrow(() -> new IllegalArgumentException("no valid input line found")).split("_").length - 1;
        List<String> outputLines = Files.lines(Paths.get(StudyConstants.getResourceAsFilepath(FILENAME)))
                .filter(ComposeItemSetMinerDataSet::filterFunctionalProteins)
                .map(ComposeItemSetMinerDataSet::handleLine)
                .filter(line -> !line.isEmpty())
                .collect(Collectors.toList());

        String output = outputLines.stream()
                .collect(Collectors.joining(System.lineSeparator(), "@RELATION ism" + System.lineSeparator() +
                        "@ATTRIBUTE id             STRING" + System.lineSeparator() +
                        IntStream.range(0, numberOfGroups)
                                .mapToObj(ComposeItemSetMinerDataSet::mapToHeaderLine)
                                .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                        "@ATTRIBUTE class          {functional,non-functional}" + System.lineSeparator() +
                        "@DATA" + System.lineSeparator(), ""));

        System.out.println(output);
    }

    private static List<String> functionalPdbIds;

    static {
        try {
            functionalPdbIds = Files.lines(Paths.get(StudyConstants.getResourceAsFilepath("studies/gmlvq/positives.txt")))
                    .map(line -> line.split("_")[0])
                    .distinct()
                    .collect(Collectors.toList());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static boolean filterFunctionalProteins(String line) {
        if(FUNCTIONAL) {
            return true;
        }

        String pdbId = line.split("_")[0];
        if(functionalPdbIds.contains(pdbId)) {
            System.out.println("filtered " + pdbId);
            return false;
        }

        return true;
    }

    private static String mapToHeaderLine(int index) {
        index++;
        return "@ATTRIBUTE energy" + index + "        NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE loopFraction" + index + "  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE rasa" + index + "          NUMERIC";
    }

    private static String handleLine(String line) {
        System.out.println(line);
        try {
            String[] sectionSplit = line.split("_");
            String pdbId = sectionSplit[0];
            Protein protein = ProteinParser.source(pdbId).parse();

            // compute features
            featureProviders.forEach(featureProvider -> featureProvider.process(protein));

            // extract peculiar residues
            List<Group> groups = Stream.of(sectionSplit)
                    // skip pdbId
                    .skip(1)
                    .map(residueSection -> extractResidue(protein, residueSection))
                    .collect(Collectors.toList());

            return groups.stream()
                    .map(ComposeItemSetMinerDataSet::mapToString)
                    .collect(Collectors.joining(",", line + ",", "," + (FUNCTIONAL ? "" : "non-") + "functional"));
        } catch (NoSuchElementException | NullPointerException | NumberFormatException e) {
            // thrown upon missing backbone atoms during secondary structure assignment
            return "";
        }
    }

    private static String mapToString(Group group) {
        return DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(EnergyProfile.class)) + "," +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(LoopFraction.class)) + "," +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(AccessibleSurfaceArea.class));
    }

    private static Group extractResidue(Protein protein, String residueSection) {
        String[] split = residueSection.split("-");
        Group group = protein.select().chainName(split[0]).residueNumber(Integer.valueOf(split[1].substring(1))).asGroup();
        // check for integrity
        if(!group.getGroupInformation().getOneLetterCode().equals(split[1].substring(0, 1))) {
            throw new IllegalArgumentException("amino acid does not match expectation: " + residueSection + " found " + group.getGroupInformation().getOneLetterCode());
        }
        return group;
    }
}

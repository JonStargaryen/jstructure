package studies.gmlvq;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * The abstract composer for GMLVQ_MAIN-data sets dealing with itemset miner data.
 * Created by bittrich on 5/22/17.
 */
public abstract class AbstractItemSetMinerDataSetComposer {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####");
    private static final List<AbstractFeatureProvider> featureProviders = Stream.of(AccessibleSurfaceArea.class,
            LoopFraction.class,
            EnergyProfile.class)
            .map(FeatureProviderRegistry::resolve)
            .collect(Collectors.toList());

    private List<String> functionalPdbIds;

    protected AbstractItemSetMinerDataSetComposer(String positiveFilename, String negativeFilename) throws IOException {
        this.functionalPdbIds = Files.lines(Paths.get(positiveFilename))
                .map(line -> line.split("_")[0])
                .distinct()
                .collect(Collectors.toList());
        handleCase(positiveFilename, true);
        handleCase(negativeFilename, false);
    }

    private void handleCase(String filename, boolean functional) throws IOException {
        // the number of considered residues
        int numberOfGroups = Files.lines(Paths.get(filename))
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException("no valid input line found"))
                .split("_").length - 1;
        List<String> outputLines = Files.lines(Paths.get(filename))
                .filter(line -> filterFunctionalProteins(line, functional))
                .map(line -> handleLine(line, functional))
                .filter(line -> !line.isEmpty())
                .collect(Collectors.toList());

        String output = outputLines.stream()
                .collect(Collectors.joining(System.lineSeparator(), "@RELATION ism" + System.lineSeparator() +
                        "@ATTRIBUTE id             STRING" + System.lineSeparator() +
                        IntStream.range(0, numberOfGroups)
                                .mapToObj(this::mapToHeaderLine)
                                .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                        "@ATTRIBUTE class          {functional,non-functional}" + System.lineSeparator() +
                        "@DATA" + System.lineSeparator(), ""));

        System.out.println(output);
    }

    private boolean filterFunctionalProteins(String line, boolean functional) {
        if(functional) {
            return true;
        }

        String pdbId = line.split("_")[0];
        if(functionalPdbIds.contains(pdbId)) {
            System.out.println("filtered " + pdbId);
            return false;
        }

        return true;
    }

    private String mapToHeaderLine(int index) {
        index++;
        return "@ATTRIBUTE energy" + index + "        NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE loopFraction" + index + "  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE rasa" + index + "          NUMERIC";
    }

    private String handleLine(String line, boolean functional) {
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

            String partialOutput = groups.stream()
                    .map(this::mapToString)
                    .collect(Collectors.joining(",", line + ",", "," + (functional ? "" : "non-") + "functional"));
            System.out.println(partialOutput);
            return partialOutput;
        } catch (NoSuchElementException | NullPointerException | NumberFormatException e) {
            e.printStackTrace();
            // thrown upon missing backbone atoms during secondary structure assignment
            return "";
        }
    }

    private String mapToString(Group group) {
        return DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(LoopFraction.class).getLoopFraction()) + "," +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea());
    }

    private Group extractResidue(Protein protein, String residueSection) {
        String[] split = residueSection.split("-");
        Group group = protein.select().chainName(split[0]).residueNumber(Integer.valueOf(split[1].substring(1))).asGroup();
        // check for integrity
        if(!group.getGroupInformation().getOneLetterCode().equals(split[1].substring(0, 1))) {
            throw new IllegalArgumentException("amino acid does not match expectation: " + residueSection + " found " + group.getGroupInformation().getOneLetterCode());
        }
        return group;
    }
}

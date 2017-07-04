package studies.gmlvq.fingerprint;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Composes data sets for the finger print data.
 * Created by bittrich on 6/23/17.
 */
public class FingerPrintDataSetComposer {
    private static final Logger logger = LoggerFactory.getLogger(FingerPrintDataSetComposer.class);
    private static final String DECOY = "decoy";
    private static final String FINGERPRINT = "fingerprint";
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####");
    private static final List<AbstractFeatureProvider> featureProviders = Stream.of(AccessibleSurfaceArea.class,
            LoopFraction.class,
            EnergyProfile.class,
            PLIPInteractionContainer.class)
            .map(FeatureProviderRegistry::resolve)
            .collect(Collectors.toList());
    private ProteinIdentifier lastProteinIdentifier = null;
    private Protein lastProtein = null;
    private static final double MINIMAL_PERCENTAGE = 0.25;
    private static final int MAX_SIZE = 1000;

    /**
     *
     * @param familyPath the path to the motif directory - supposed to contain 2 subdirectories: decoy and fingerprint
     */
    public FingerPrintDataSetComposer(Path familyPath, Path outputDirectory) {
        String familyName = familyPath.toFile().getName();
        logger.info("starting job on {}", familyName);
        Path parentPath = familyPath.getParent();
        Path resultPath = parentPath.resolve("results").resolve(familyName);

        try {
            Files.list(resultPath)
                    .forEach(motifPath -> {
                        try {
                            String motifName = motifPath.toFile().getName();
                            Path outputPath = outputDirectory.resolve(familyName + "_" + motifName + ".arff");
                            if(Files.exists(outputPath)) {
                                logger.info("skipping already present results for {}/{}", familyName, motifName);
                                return;
                            }

                            logger.info("starting sub-job on {}/{}", familyName, motifName);
                            Path decoyDirectory = motifPath.resolve(DECOY);
                            Path fingerprintDirectory = motifPath.resolve(FINGERPRINT);

                            double decoyCount = Files.list(decoyDirectory).count();
                            double fingerprintCount = Files.list(fingerprintDirectory).count();
                            double nonRedundantCount = Files.list(familyPath.resolve("nr")).count();

                            if(fingerprintCount < nonRedundantCount * MINIMAL_PERCENTAGE) {
                                // skip sparsely populated results
                                logger.warn("skipping sparsely populated results for {}/{}", familyName, motifName);
                                return;
                            }

                            if(decoyCount > MAX_SIZE) {
                                logger.warn("encountering highly populated decoy set for {}/{} - sampling down to {}", familyName, motifName, MAX_SIZE);
                            }

                            int numberOfGroups = Files.list(decoyDirectory)
                                    .findFirst()
                                    .get()
                                    .toFile()
                                    .getName()
                                    .split("_")
                                    .length - 1;
                            String headerString = "@RELATION ism" + System.lineSeparator() +
                                    "@ATTRIBUTE id             STRING" + System.lineSeparator() +
                                    IntStream.range(0, numberOfGroups)
                                            .mapToObj(this::mapToHeaderLine)
                                            .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                                    "@ATTRIBUTE class          {functional,non-functional}" + System.lineSeparator() +
                                    "@DATA";
                            List<Path> decoyPaths = Files.list(decoyDirectory)
                                    // max size of decoy set
                                    .limit(MAX_SIZE)
                                    .sorted(Path::compareTo)
                                    .collect(Collectors.toList());
                            String decoyString = decoyPaths.stream()
                                    .map(this::handleLine)
                                    .filter(Optional::isPresent)
                                    .map(Optional::get)
                                    .peek(logger::info)
                                    .collect(Collectors.joining(System.lineSeparator()));
                            String fingerprintString = Files.list(fingerprintDirectory)
                                    .map(this::handleLine)
                                    .filter(Optional::isPresent)
                                    .map(Optional::get)
                                    .peek(logger::info)
                                    .collect(Collectors.joining(System.lineSeparator()));

                            Files.write(outputPath,
                                    (headerString + System.lineSeparator() + decoyString + System.lineSeparator() + fingerprintString).getBytes());
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    });
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private String mapToHeaderLine(int index) {
        index++;
        return "@ATTRIBUTE energy" + index + "        NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE loopFraction" + index + "  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE rasa" + index + "          NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE halogen" + index + "       NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE hydrogen" + index + "      NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE metalComplex" + index + "  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE piCation" + index + "      NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE piStackings" + index + "   NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE saltBridges" + index + "   NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE waterBridges" + index + "  NUMERIC";
    }

    private Optional<String> handleLine(Path path) {
        try {
            String filename = path.toFile().getName().split("\\.")[0];

            String[] sectionSplit = filename.split("_");
            ProteinIdentifier proteinIdentifier = ProteinIdentifier.createFromPdbId(filename.split("_")[0]);
            Protein protein;
            if(proteinIdentifier.equals(lastProteinIdentifier)) {
                protein = lastProtein;
                if(lastProtein == null) {
                    return Optional.empty();
                }
            } else {
                lastProteinIdentifier = proteinIdentifier;
                try {
                    protein = ProteinParser.source(proteinIdentifier.getPdbId()).parse();
                    featureProviders.forEach(featureProvider -> featureProvider.process(protein));
                    lastProtein = protein;
                } catch (Exception e) {
                    lastProtein = null;
                    return Optional.empty();
                }
            }

            // extract peculiar residues
            List<Group> groups = Stream.of(sectionSplit)
                    // skip pdbId
                    .skip(1)
                    .map(residueSection -> extractResidue(protein, residueSection))
                    .collect(Collectors.toList());

            PLIPInteractionContainer container = protein.getFeatureContainer().getFeature(PLIPInteractionContainer.class);

            // determine if functional or non-functional record
            boolean functional = path.getParent().toFile().getName().equals("fingerprint");

            return Optional.of(groups.stream()
                    .map(group -> mapToString(container, group))
                    .collect(Collectors.joining(",", filename + ",", "," + (functional ? "" : "non-") + "functional")));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private String mapToString(PLIPInteractionContainer container, Group group) {
        return DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(LoopFraction.class).getLoopFraction()) + "," +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                composeInteractionString(container, group);
    }

    private Group extractResidue(Protein protein, String residueSection) {
        String[] split = residueSection.split("-");
        Group group = protein
                .select()
                .chainName(split[0])
                .residueNumber(Integer.valueOf(split[2]))
                .asGroup();
        // check for integrity
        if(!group.getGroupPrototype().getThreeLetterCode().equalsIgnoreCase(split[1])) {
            // happens for alternative positions
            throw new IllegalArgumentException("amino acid does not match expectation: " + residueSection + " found " +
                    group.getGroupPrototype().getOneLetterCode());
        }
        return group;
    }

    private String composeInteractionString(PLIPInteractionContainer plipInteractionContainer, Group group) {
        long halogenBonds = plipInteractionContainer.getHalogenBonds().stream()
                .filter(halogenBond -> describesGroup(halogenBond, group))
                .count();
        long hydrogenBonds = plipInteractionContainer.getHydrogenBonds().stream()
                .filter(hydrogenBond -> describesGroup(hydrogenBond, group))
                .count();
        long metalComplexes = plipInteractionContainer.getMetalComplexes().stream()
                .filter(metalComplex -> describesGroup(metalComplex, group))
                .count();
        long piCation = plipInteractionContainer.getPiCationInteractions().stream()
                .filter(piCationInteraction -> describesGroup(piCationInteraction, group))
                .count();
        long piStacking = plipInteractionContainer.getPiStackings().stream()
                .filter(piStackingInteraction -> describesGroup(piStackingInteraction, group))
                .count();
        long saltBridges = plipInteractionContainer.getSaltBridges().stream()
                .filter(saltBridge -> describesGroup(saltBridge, group))
                .count();
        long waterBridges = plipInteractionContainer.getWaterBridges().stream()
                .filter(waterBridge -> describesGroup(waterBridge, group))
                .count();
        return halogenBonds + "," +
                hydrogenBonds + "," +
                metalComplexes + "," +
                piCation + "," +
                piStacking + "," +
                saltBridges + "," +
                waterBridges;
    }

    private boolean describesGroup(PLIPInteraction plipInteraction, Group group) {
        return plipInteraction.getPartner1().equals(group) || plipInteraction.getPartner2().equals(group);
    }
}

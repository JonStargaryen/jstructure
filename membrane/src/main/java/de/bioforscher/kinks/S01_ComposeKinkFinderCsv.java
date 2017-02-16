package de.bioforscher.kinks;

import de.bioforscher.Constants;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.interaction.HydrogenBond;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;

import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Compose all results into one csv.
 * Created by bittrich on 2/14/17.
 */
@Deprecated
public class S01_ComposeKinkFinderCsv {
    private static final AbstractFeatureProvider asaCalculator = new AccessibleSurfaceAreaCalculator();
    private static final AbstractFeatureProvider plipAnnotator = new PLIPAnnotator();
    private static final AbstractFeatureProvider sequenceMotifAnnotator = new SequenceMotifAnnotator();
    static final String DELIMITER = "\t";
    static final List<String> ATTRIBUTE_NAMES = Stream.of("pdbId", "chainId", "start" , "end", "kink", "gly", "pro",
            "rasaL4", "rasaL3", "rasaL2", "rasaL1", "rasa", "rasaR1","rasaR2", "rasaR3", "rasaR4", "plipL4", "plipL3",
            "plipL2", "plipL1", "plip", "plipR1", "plipR2", "plipR3", "plipR4",  "hbondL4", "hbondL3", "hbondL2",
            "hbondL1", "hbond", "hbondR1", "hbondR2", "hbondR3", "hbondR4", "sideL4", "sideL3", "sideL2", "sideL1",
            "side", "sideR1", "sideR2", "sideR3", "sideR4", "motifs", "error", "angle").collect(Collectors.toList());

    public static void main(String[] args) {
        String header = ATTRIBUTE_NAMES.stream().collect(Collectors.joining(DELIMITER));

        String output = Constants.lines(Paths.get(Constants.PDBTM_ALPHA_NR_LIST))
                .filter(Constants.isCommentLine.negate())
                .map(line -> line.substring(0, 4))
                .distinct()
                // ensure KinkFinder results exist
                .filter(pdbId -> Constants.list(Paths.get(Constants.KINK_FINDER_RESULT_PATH))
                        .anyMatch(path -> path.toFile().getName().startsWith(pdbId)))
                .map(S01_ComposeKinkFinderCsv::parseProtein)
                // when plip annotation is missing null is returned, skip those cases
                .filter(Objects::nonNull)
                .map(S01_ComposeKinkFinderCsv::composeCsvLines)
                .collect(Collectors.joining(System.lineSeparator(), header + System.lineSeparator(), ""));

        Constants.write(Paths.get(Constants.STATISTICS_PATH + "kink_finder_results_all.csv"), output.getBytes());
    }

    private static String composeCsvLines(Protein protein) {
        List<KinkFinderHelix> helices = protein.getFeatureAsList(KinkFinderHelix.class, KinkFinderParser.KINK_FINDER_ANNOTATION);
        return helices.stream()
                // filter for non-redundant
                .filter(Objects::nonNull)
                .map(kinkFinderHelix -> {
                    try {
                        return composeCsvLine(kinkFinderHelix, protein);
                    } catch (NullPointerException e) {
                        return null;
                    }
                })
                .filter(Objects::nonNull)
                .peek(System.out::println)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static String composeCsvLine(KinkFinderHelix kinkFinderHelix, Protein protein) {
        Map<String, Long> occurrences = Arrays.stream(kinkFinderHelix.getSequence().split(""))
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));

        String pdbId = kinkFinderHelix.getPdbCode().substring(0, 4);
        String chainId = kinkFinderHelix.getPdbCode().substring(4);
        int kinkPosition = kinkFinderHelix.getKinkPosition();
        int helixStart = kinkFinderHelix.getHelixStart();
        int helixEnd = kinkFinderHelix.getHelixEnd();
        long glycineCount;
        try {
            glycineCount = occurrences.get("G");
        } catch (NullPointerException e) {
            glycineCount = 0;
        }
        long prolineCount;
        try {
            prolineCount = occurrences.get("P");
        } catch (NullPointerException e) {
            prolineCount = 0;
        }
        double estimatedError = kinkFinderHelix.getEstimatedError();
        double kinkAngle = kinkFinderHelix.getKinkAngle();
        GroupContainer groupContainer = Selection.on(protein)
                .chainName(chainId)
                .residueNumber(new IntegerRange(kinkPosition - 4, kinkPosition + 4))
                .asGroupContainer();

        String sequenceMotifs;

        try {
            sequenceMotifs = Selection.on(protein)
                    .chainName(chainId)
                    .residueNumber(kinkPosition)
                    .asGroup()
                    .getFeatureAsList(SequenceMotif.class, SequenceMotifAnnotator.SEQUENCE_MOTIF)
                    .stream()
                    .map(SequenceMotif::getMotifDefinition)
                    .map(SequenceMotifDefinition::name)
                    .collect(Collectors.joining(",", "[", "]"));
        } catch (NullPointerException e) {
            sequenceMotifs = "[]";
        }

        return pdbId + DELIMITER + chainId + DELIMITER + helixStart + DELIMITER + helixEnd + DELIMITER + kinkPosition + DELIMITER + glycineCount +
                DELIMITER + prolineCount + DELIMITER +
                groupContainer.groups()
                        .mapToDouble(group -> group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA))
                        .mapToObj(Constants.DECIMAL_FORMAT::format)
                        .collect(Collectors.joining(DELIMITER))
                + DELIMITER +
                groupContainer.groups()
                        .map(group -> group.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS))
                        .mapToInt(Collection::size)
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(DELIMITER))
                + DELIMITER +
                groupContainer.groups()
                        .map(group -> group.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS)
                                .stream()
                                .filter(plipInteraction -> plipInteraction instanceof HydrogenBond)
                                // ensure that interactions describe backbone interactions TODO decide actually on annotated, interacting atoms
                                .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) < 7)
                                .count()
                        )
                        .map(String::valueOf)
                        .collect(Collectors.joining(DELIMITER))
                + DELIMITER +
                groupContainer.groups()
                        .map(group -> group.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS)
                                .stream()
                                // ensure that interactions describe backbone interactions
                                .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) > 7)
                                .count()
                        )
                        .map(String::valueOf)
                        .collect(Collectors.joining(DELIMITER))
                + DELIMITER + sequenceMotifs + DELIMITER + Constants.DECIMAL_FORMAT.format(estimatedError) + DELIMITER + Constants.DECIMAL_FORMAT.format(kinkAngle);
    }

    private static Protein parseProtein(String pdbId) {
        Protein protein = ProteinParser.parsePDBFile(Constants.STRUCTURE_PATH + pdbId + Constants.PDB_SUFFIX);

        KinkFinderParser.parseKinkFinderFile(protein, Paths.get(Constants.KINK_FINDER_RESULT_PATH + pdbId + ".kinks"));
        // assign additional features
        plipAnnotator.process(protein);
        asaCalculator.process(protein);
        sequenceMotifAnnotator.process(protein);

        if(protein.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS).size() == 0) {
            return null;
        }

        return protein;
    }
}

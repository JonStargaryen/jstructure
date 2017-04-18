package bioforscher.kinks;

import bioforscher.Constants;
import bioforscher.opm.OPMParser;
import bioforscher.opm.TMHelix;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;

import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 *
 * Created by bittrich on 2/16/17.
 */
public class S01_ComposeKinkFinderData {
    private static final AbstractFeatureProvider asaCalculator = new AccessibleSurfaceAreaCalculator();
    private static final AbstractFeatureProvider plipAnnotator = new PLIPAnnotator();
    private static final AbstractFeatureProvider sequenceMotifAnnotator = new SequenceMotifAnnotator();
    private static final String DELIMITER = "\t";

    public static void main(String[] args) {
        List<String> nrpdbIds = Constants.lines(Paths.get(Constants.PDBTM_ALPHA_NR_LIST))
                .filter(Constants.isCommentLine.negate())
                .skip(2)
                // ensure KinkFinder results exist
                .filter(line -> Constants.list(Paths.get(Constants.KINK_FINDER_RESULT_PATH))
                        .anyMatch(path -> path.toFile().getName().startsWith(line.substring(0, 4))))
                // ensure OPM annotation exists
                .filter(line -> Constants.list(Paths.get(Constants.OPM_PATH))
                        .anyMatch(path -> path.toFile().getName().startsWith(line.substring(0, 4))))
                .limit(5)
                .collect(Collectors.toList());

        // raw list of proteins, but not necessarily non-redundant chains
        List<Protein> proteins = nrpdbIds.stream()
                .map(line -> line.substring(0, 4))
                .map(pdbId -> Constants.STRUCTURE_PATH + pdbId + Constants.PDB_SUFFIX)
                .distinct()
                .map(Paths::get)
                .map(path -> ProteinParser.source(path).parse())
                .peek(plipAnnotator::process)
                .peek(sequenceMotifAnnotator::process)
                .collect(Collectors.toList());

        StringJoiner output = new StringJoiner(System.lineSeparator());
        nrpdbIds.forEach(nrpdbId -> {
            String pdbId = nrpdbId.substring(0, 4);
            String chainId = nrpdbId.substring(5);
            Protein protein = proteins.stream()
                    .filter(p -> p.getName().equalsIgnoreCase(pdbId))
                    .findFirst()
                    .get();
            Chain chain = Selection.on(protein)
                    .chainName(chainId)
                    .asChain();

            // ensure PLIP annotation is present for the wanted chain
            PLIPInteractionContainer plipInteractionContainer = protein.getFeature(PLIPInteractionContainer.class, PLIPAnnotator.PLIP_INTERACTIONS);
            if(plipInteractionContainer.getInteractions().isEmpty()) {
                return;
            }

            List<PLIPInteraction> interactionsInChain = plipInteractionContainer.getInteractions().stream()
                    .filter(plipInteraction -> plipInteraction.getPartner1().getParentChain().getChainId().equals(chainId))
                    .collect(Collectors.toList());
            if(interactionsInChain.isEmpty()) {
                return;
            }

            OPMParser.parse(protein, Paths.get(Constants.OPM_PATH + pdbId + ".opm_"));
            List<TMHelix> helices = protein.getFeatureAsList(TMHelix.class, OPMParser.TM_HELIX);

            KinkFinderParser.parseKinkFinderFile(protein, Paths.get(Constants.KINK_FINDER_RESULT_PATH + pdbId + ".kinks"));
            List<KinkFinderHelix> kinks = new ArrayList<>(protein.getFeatureAsList(KinkFinderHelix.class, KinkFinderParser.KINK_FINDER_ANNOTATION));
            List<KinkFinderHelix> significantKinks = kinks.stream()
                    .filter(kinkFinderHelix -> kinkFinderHelix.getKinkAngle() > 20.0)
                    .collect(Collectors.toList());
            int kinkCount = significantKinks.size();

            System.out.println(pdbId + "-" + chainId + ": " + helices.size() + " helices, " + interactionsInChain.size() + " interactions, " + kinkCount + " kinks");
            asaCalculator.process(protein);

            List<KinkFinderHelix> nonKinks = new ArrayList<>();
            // for the given structure select as many non-kink observations
            for (KinkFinderHelix referenceKink : significantKinks) {
                try {
                    int nonKinkPosition = generateRandomNonKinkPosition(helices, significantKinks, nonKinks);
                    int[] helixPositions = determineHelix(nonKinkPosition, helices);
                    String sequence = Selection.on(chain)
                            .residueNumber(new IntegerRange(nonKinkPosition - 6, nonKinkPosition + 6))
                            .asGroupContainer()
                            .getAminoAcidSequence();

                    nonKinks.add(new KinkFinderHelix(referenceKink.getPdbCode(),
                            helixPositions[0],
                            helixPositions[1],
                            nonKinkPosition - 6,
                            nonKinkPosition + 6,
                            nonKinkPosition,
                            0.0,
                            sequence,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0));
                } catch (ArrayIndexOutOfBoundsException e) {
                    System.out.println("no more positions available for " + pdbId);
                } catch (NoSuchElementException e) {
                    System.out.println(e.getLocalizedMessage());
                }
            }

            Stream.of(significantKinks, nonKinks)
                    .flatMap(Collection::stream)
                    .map(kinkFinderHelix -> {
                        Map<String, Long> occurrences = Arrays.stream(kinkFinderHelix.getSequence().split(""))
                                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));

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
//                                + DELIMITER +
//                                groupContainer.groups()
//                                        .map(group -> group.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS))
//                                        .mapToInt(Collection::size)
//                                        .mapToObj(String::valueOf)
//                                        .collect(Collectors.joining(DELIMITER))
//                                + DELIMITER +
//                                groupContainer.groups()
//                                        .map(group -> group.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS)
//                                                .stream()
//                                                .filter(plipInteraction -> plipInteraction instanceof HydrogenBond)
//                                                // ensure that interactions describe backbone interactions TODO decide actually on annotated, interacting atoms
//                                                .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) < 7)
//                                                .count()
//                                        )
//                                        .map(String::valueOf)
//                                        .collect(Collectors.joining(DELIMITER))
//                                + DELIMITER +
//                                groupContainer.groups()
//                                        .map(group -> group.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS)
//                                                .stream()
//                                                // ensure that interactions describe side chain interactions
//                                                .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) > 7)
//                                                .count()
//                                        )
//                                        .map(String::valueOf)
//                                        .collect(Collectors.joining(DELIMITER))
                                + DELIMITER + sequenceMotifs + DELIMITER + (kinkFinderHelix.isSignificantKink() ? "K" : "n") + DELIMITER + Constants.DECIMAL_FORMAT.format(kinkAngle);
                    })
                    .forEach(output::add);
        });

        System.out.println(output.toString());
    }

    private static int[] determineHelix(int nonKinkPosition, List<TMHelix> helices) {
        TMHelix tmHelix = helices.stream()
                .filter(h -> h.getStartGroup().getResidueNumber() <= nonKinkPosition && h.getEndGroup().getResidueNumber() >= nonKinkPosition)
                .findFirst()
                .orElseThrow(() -> new NoSuchElementException("no helix annotated around " + nonKinkPosition));

        return new int[] { tmHelix.getStartGroup().getResidueNumber(), tmHelix.getEndGroup().getResidueNumber() };
    }

    /**
     * Generate a random index which can represent a non-kink position.
     * Their should be neither another kink or a previously selected position in the direct neighborhood. It has to be
     * located in a TMHelix.
     * @param helices what positions are available?
     * @param kinks and which are taken?
     * @param nonKinks also taken?
     * @return a random index which is still available
     */
    private static int generateRandomNonKinkPosition(List<TMHelix> helices, List<KinkFinderHelix> kinks, List<KinkFinderHelix> nonKinks) {
        final int minimalSeparationToOtherKink = 2;
        final int minimalSeparationToHelixEnd = 6;

        List<Integer> availablePositions = helices.stream()
                .flatMap(tmHelix -> IntStream.range(tmHelix.getStartGroup().getResidueNumber() + minimalSeparationToHelixEnd, tmHelix.getEndGroup().getResidueNumber() - minimalSeparationToHelixEnd)
                        .boxed())
                .collect(Collectors.toList());

        List<Integer> positionsToIgnore = Stream.of(kinks, nonKinks)
                .flatMap(Collection::stream)
                .map(KinkFinderHelix::getKinkPosition)
                .flatMap(position -> IntStream.range(position - minimalSeparationToOtherKink, position + minimalSeparationToOtherKink + 1)
                        .boxed())
                .distinct()
                .collect(Collectors.toList());

        availablePositions.removeAll(positionsToIgnore);

        Collections.shuffle(availablePositions);

        return availablePositions.get(0);
    }
}

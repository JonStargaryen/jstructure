package studies.gmlvq.ahah;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import studies.StudyConstants;

import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Classify helices containing kinks based on the AHAH data set.
 * Created by bittrich on 6/27/17.
 */
public class S01_ClassifyKinksOnAhahDataSet {
    private static final Logger logger = LoggerFactory.getLogger(S01_ClassifyKinksOnAhahDataSet.class);
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####");
    private static final List<FeatureProvider> featureProviders = Stream.of(new AccessibleSurfaceAreaCalculator(),
            new EnergyProfileCalculator(),
            new PLIPIntraMolecularAnnotator())
            .collect(Collectors.toList());
    private static Map<String, Optional<Structure>> proteinMap;

    public static void main(String[] args) {
        proteinMap = StudyConstants.lines(StudyConstants.AHAH_DATA_SET)
                .filter(line -> !line.startsWith(">") && !line.startsWith("helix_id"))
                .map(line -> line.split("_")[0])
                .distinct()
                .collect(Collectors.toMap(Function.identity(),
                        S01_ClassifyKinksOnAhahDataSet::handlePdbId));

        String output = StudyConstants.lines(StudyConstants.AHAH_DATA_SET)
                .filter(line -> !line.startsWith(">") && !line.startsWith("helix_id"))
                .map(S01_ClassifyKinksOnAhahDataSet::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .peek(System.out::println)
                .collect(Collectors.joining(System.lineSeparator(),
                        "@RELATION ism" + System.lineSeparator() +
                                "@ATTRIBUTE id            STRING" + System.lineSeparator() +
                                "@ATTRIBUTE energy        NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE rasa          NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE globalInt     NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE helixHydrogen NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE longRangeInt  NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE backboneInt   NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE sideChainInt  NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE mixedInt      NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE glycine       NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE proline       NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE missingHydro3 NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE missingHydro4 NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE missingHydro5 NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE class         {KINK,CURVED,STRAIGHT,UNCLASSIFIED}" + System.lineSeparator() +
                                "@ATTRIBUTE kink          {KINK,NONE}" + System.lineSeparator() +
                                "@DATA" + System.lineSeparator(),
                        System.lineSeparator()));

        StudyConstants.write(StudyConstants.AHAH_DATA_SET.getParent().resolve("ahah.arff"), output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            String[] commaSplit = line.split(",");
            String[] idSplit = commaSplit[0].split("_");

            Structure protein = proteinMap.get(idSplit[0]).get();
            Chain chain = protein.select()
                    .chainName(idSplit[1])
                    .asChain();
            IntegerRange range = new IntegerRange(Integer.valueOf(idSplit[2]), Integer.valueOf(idSplit[3]));
            GroupContainer groups = chain.select()
                    .aminoAcids()
                    .residueNumber(range)
                    .asIsolatedStructure();
            double size = groups.aminoAcids().count();

            double averageEnergy = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(EnergyProfile.class))
                    .mapToDouble(EnergyProfile::getSolvationEnergy)
                    .average()
                    .getAsDouble();
            double averageRasa = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(AccessibleSurfaceArea.class))
                    .mapToDouble(AccessibleSurfaceArea::getRelativeAccessibleSurfaceArea)
                    .average()
                    .getAsDouble();
            // all hydrogen bonds observable for residues of this helix
            double globalInteractions = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getInteractions)
                    .mapToInt(List::size)
                    .sum() / size;
            // all hydrogen bonds occurring within the helix exclusively
            double helixHydrogenBonds =  groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getHydrogenBonds)
                    .flatMap(Collection::stream)
                    .filter(hydrogenBond -> groups.getGroups().contains(hydrogenBond.getPartner1()) && groups.getGroups().contains(hydrogenBond.getPartner2()))
                    .count() / size;
            // all interactions which are only partly within the helix
            double longRangeInteractions = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getInteractions)
                    .flatMap(Collection::stream)
                    .filter(interaction -> (groups.getGroups().contains(interaction.getPartner1()) && !groups.getGroups().contains(interaction.getPartner2())) || (!groups.getGroups().contains(interaction.getPartner1()) && groups.getGroups().contains(interaction.getPartner2())))
                    .count() / size;
            // of global interactions: how many are backbone interactions
            double backboneInteractions = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getInteractions)
                    .flatMap(Collection::stream)
                    .filter(PLIPInteraction::isBackboneInteraction)
                    .count() / size;
            // of global interactions: how many are side-chain interactions
            double sideChainInteractions = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getInteractions)
                    .flatMap(Collection::stream)
                    .filter(PLIPInteraction::isSideChainInteraction)
                    .count() / size;
            // of global interactions: how many are mixed interactions
            double mixedInteractions = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getInteractions)
                    .flatMap(Collection::stream)
                    .filter(PLIPInteraction::isMixedInteraction)
                    .count() / size;

            double glycines = groups.getAminoAcidSequence().split("G").length - 1 / size;
            double prolines = groups.getAminoAcidSequence().split("P").length - 1 / size;

            // count missing backbone interactions
            double missingBackboneInteractions3 = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getHydrogenBonds)
                    .flatMap(Collection::stream)
                    // first interacting atom should be 'O'
                    .filter(interation -> "O".equals(interation.getAcceptor().getName()))
                    // donor should be backbone nitrogen
                    .filter(interaction -> "N".equals(interaction.getDonor().getName()))
                    // should be i -> i + 3
                    .filter(interaction -> interaction.getPartner1().getResidueIdentifier().getResidueNumber() + 3 == interaction.getPartner2().getResidueIdentifier().getResidueNumber())
                    .count() / size;
            double missingBackboneInteractions4 = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getHydrogenBonds)
                    .flatMap(Collection::stream)
                    // first interacting atom should be 'O'
                    .filter(interation -> "O".equals(interation.getAcceptor().getName()))
                    // donor should be backbone nitrogen
                    .filter(interaction -> "N".equals(interaction.getDonor().getName()))
                    // should be i -> i + 4
                    .filter(interaction -> interaction.getPartner1().getResidueIdentifier().getResidueNumber() + 4 == interaction.getPartner2().getResidueIdentifier().getResidueNumber())
                    .count() / size;
            double missingBackboneInteractions5 = groups.aminoAcids()
                    .map(AminoAcid::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getHydrogenBonds)
                    .flatMap(Collection::stream)
                    // first interacting atom should be 'O'
                    .filter(interation -> "O".equals(interation.getAcceptor().getName()))
                    // donor should be backbone nitrogen
                    .filter(interaction -> "N".equals(interaction.getDonor().getName()))
                    // should be i -> i + 5
                    .filter(interaction -> interaction.getPartner1().getResidueIdentifier().getResidueNumber() + 5 == interaction.getPartner2().getResidueIdentifier().getResidueNumber())
                    .count() / size;

            return Optional.of(
                    // id
                    commaSplit[0] + "," +
                    DECIMAL_FORMAT.format(averageEnergy) + "," +
                    DECIMAL_FORMAT.format(averageRasa) + "," +
                    DECIMAL_FORMAT.format(globalInteractions) + "," +
                    DECIMAL_FORMAT.format(helixHydrogenBonds) + "," +
                    DECIMAL_FORMAT.format(longRangeInteractions) + "," +
                    DECIMAL_FORMAT.format(backboneInteractions) + "," +
                    DECIMAL_FORMAT.format(sideChainInteractions) + "," +
                    DECIMAL_FORMAT.format(mixedInteractions) + "," +
                    DECIMAL_FORMAT.format(glycines) + "," +
                    DECIMAL_FORMAT.format(prolines) + "," +
                    DECIMAL_FORMAT.format(missingBackboneInteractions3) + "," +
                    DECIMAL_FORMAT.format(missingBackboneInteractions4) + "," +
                    DECIMAL_FORMAT.format(missingBackboneInteractions5) + "," +
                    // 4 state AHAH class
                    commaSplit[1] + "," +
                    (commaSplit[1].equals("KINK") ? "KINK" : "NONE")
            );
        } catch (Exception e) {
            return Optional.empty();
        }
    }

    private static Optional<Structure> handlePdbId(String pdbId) {
        try {
            logger.info("fetching and annotating {}", pdbId);
            Structure protein = StructureParser.source(pdbId).parse();
            featureProviders.forEach(featureProvider -> featureProvider.process(protein));
            return Optional.of(protein);
        } catch (Exception e) {
            return Optional.empty();
        }
    }
}

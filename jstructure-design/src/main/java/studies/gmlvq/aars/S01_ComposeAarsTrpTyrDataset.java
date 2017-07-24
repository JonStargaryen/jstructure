package studies.gmlvq.aars;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPLigandAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static de.bioforscher.jstructure.StandardFormat.format;

public class S01_ComposeAarsTrpTyrDataset {
    private static final AccessibleSurfaceAreaCalculator accessibleSurfaceAreaCalculator = new AccessibleSurfaceAreaCalculator();
    private static final EnergyProfileCalculator energyProfileCalculator = new EnergyProfileCalculator();
    private static final LoopFractionCalculator loopFractionCalculator = new LoopFractionCalculator();
    private static final PLIPLigandAnnotator plipAnnotator = new PLIPLigandAnnotator();

    /**
     * All residues with interactions according to PLIP minus that which are annotated in the binding mode matrix.
     */
    public enum ResidueNumbers {
        TRP(269,273,276,287,311,313,316,317,477,1200,1204,1207,1208,1295,1306,1416,1417,1436,1440),
        TYR(269,270,272,311,313,316,321,337,343,398,401,406,1200,1204,1207,1417);

        private int[] residueNumbers;

        ResidueNumbers(int... residueNumbers) {
            this.residueNumbers = residueNumbers;
        }

        public int[] getResidueNumbers() {
            return residueNumbers;
        }
    }

    public static void main(String[] args) {
        Map<String, List<String>> structureMap = TestUtils.getResourceAsStream("aars/ids.list")
                .filter(line -> !line.startsWith("id"))
                .collect(Collectors.groupingBy(line -> line.split("_")[0]));

        String output = structureMap.entrySet().stream()
                .flatMap(entry -> handleBin(entry.getKey(), entry.getValue()))
                .collect(Collectors.joining(System.lineSeparator(),
                        "@RELATION trptyr" + System.lineSeparator() +
                                "@ATTRIBUTE id              STRING" + System.lineSeparator() +
                                "@ATTRIBUTE representative  STRING" + System.lineSeparator() +
                                "@ATTRIBUTE bindingSiteSize NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE rasa            NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE energy          NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE loopfraction    NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE halogenBonds    NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE hydrogenBonds   NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE metalComplexes  NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE piCation        NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE piStacking      NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE saltBridges     NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE waterBridges    NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE amide           NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE amino           NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE carboxylate     NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE guanidinium     NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE hydroxyl        NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE imidizole       NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE none            NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE thiol           NUMERIC" + System.lineSeparator() +
                                "@ATTRIBUTE class           {TRP,TYR}" + System.lineSeparator() +
                                "@DATA" + System.lineSeparator(),
                        ""));

        System.out.println(output);
    }

    private static Stream<String> handleBin(String pdbId, List<String> lines) {
        Structure structure = StructureParser.source(pdbId).parse();
        accessibleSurfaceAreaCalculator.process(structure);
        energyProfileCalculator.process(structure);
        loopFractionCalculator.process(structure);
        plipAnnotator.process(structure);

        return lines.stream()
                .map(line -> handleLine(structure, line))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .peek(System.out::println);
    }

    private static Optional<String> handleLine(Structure structure, String line) {
        try {
            String[] split = line.split(",");
            ResidueNumbers type = ResidueNumbers.valueOf(split[1].toUpperCase());

            Chain chain = structure.select()
                    .chainName(split[0].split("_")[1])
                    .asChain();

            List<AminoAcid> selectedAminoAcids = getSelectedResidues(chain, type.getResidueNumbers());
            double bindingSiteSize = selectedAminoAcids.size();

            double rasa = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(AccessibleSurfaceArea.class))
                    .mapToDouble(AccessibleSurfaceArea::getRelativeAccessibleSurfaceArea)
                    .average()
                    .getAsDouble();
            double energy = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(EnergyProfile.class))
                    .mapToDouble(EnergyProfile::getSolvationEnergy)
                    .average()
                    .getAsDouble();
            double loopFraction = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(LoopFraction.class))
                    .mapToDouble(LoopFraction::getLoopFraction)
                    .average()
                    .getAsDouble();

            double halogenBonds = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getHalogenBonds)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;
            double hydrogenBonds = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getHydrogenBonds)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;
            double metalComplexes = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getMetalComplexes)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;
            double piCationInteractions = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getPiCationInteractions)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;
            double piStackings = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getPiStackings)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;
            double saltBridges = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getSaltBridges)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;
            double waterBridges = selectedAminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getWaterBridges)
                    .mapToInt(Collection::size)
                    .sum() / bindingSiteSize;

            Map<GroupPrototype.GutteridgeGrouping, Long> gutteridgeCounts = selectedAminoAcids.stream()
                    .map(AminoAcid::getGroupPrototype)
                    .collect(Collectors.groupingBy(GroupPrototype::getGutteridgeGrouping, Collectors.counting()));

            return Optional.of(split[0] + "," +
                    split[2] + "," +
                    format(bindingSiteSize) + "," +
                    format(rasa) + "," +
                    format(energy) + "," +
                    format(loopFraction) + "," +
                    format(halogenBonds) + "," +
                    format(hydrogenBonds) + "," +
                    format(metalComplexes) + "," +
                    format(piCationInteractions) + "," +
                    format(piStackings) + "," +
                    format(saltBridges)  + "," +
                    format(waterBridges) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.AMIDE, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.AMINO, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.CARBOXYLATE, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.GUANIDINIUM, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.HYDROXYL, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.IMIDAZOLE, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.NONE, 0L) / bindingSiteSize) + "," +
                    format(gutteridgeCounts.getOrDefault(GroupPrototype.GutteridgeGrouping.THIOL, 0L) / bindingSiteSize) + "," +
                    type.name());
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static List<AminoAcid> getSelectedResidues(Chain chain, int[] residueNumbers) {
        return chain.select()
                .aminoAcids()
                .residueNumber(residueNumbers)
                .asFilteredGroups()
                .map(AminoAcid.class::cast)
                .collect(Collectors.toList());
    }
}

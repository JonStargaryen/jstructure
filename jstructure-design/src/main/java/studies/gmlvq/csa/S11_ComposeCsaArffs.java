package studies.gmlvq.csa;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.mapping.SiftsMappingAnnotator;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import studies.StudyConstants;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class S11_ComposeCsaArffs {
    private static final AccessibleSurfaceAreaCalculator ACCESSIBLE_SURFACE_AREA_CALCULATOR = new AccessibleSurfaceAreaCalculator();
    private static final DictionaryOfProteinSecondaryStructure DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE = new DictionaryOfProteinSecondaryStructure();
    private static final EnergyProfileCalculator ENERGY_PROFILE_CALCULATOR = new EnergyProfileCalculator();
    private static final LoopFractionCalculator LOOP_FRACTION_CALCULATOR = new LoopFractionCalculator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) throws IOException {
        Files.list(Paths.get("/home/bittrich/git/gmlvq_main/data/csa/csa_matches/"))
                .forEach(S11_ComposeCsaArffs::handleDirectory);
    }

    private static void handleDirectory(Path csaMatchPath) {
        String motifName = csaMatchPath.toFile().getName();
        Path outputPath = csaMatchPath.getParent()
                .getParent()
                .resolve("arff")
                .resolve(motifName + ".arff");

        String pdbId = motifName.split("_")[0];
        String chainId = motifName.split("_")[1].split("-")[0];
        String ecNumber = SiftsMappingAnnotator.getEnzymeLineForPdbId(pdbId, chainId).get()[3];

        System.out.println(pdbId + "_" + chainId + " : " + ecNumber);

        BufferedWriter bufferedWriter = StudyConstants.newBufferedWriter(outputPath);
        StudyConstants.list(csaMatchPath)
                .map(match -> handleMatch(match, ecNumber))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .peek(System.out::println)
                .forEach(line -> {
                    try {
                        bufferedWriter.write(line);
                        bufferedWriter.newLine();
                        bufferedWriter.flush();
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }

    private static Optional<String> handleMatch(Path path, String originalEcNumber) {
        String motifId = path.toFile().getName().split(".pdb")[0];
        String[] filenameSplit = motifId.split("_");
        String pdbId = filenameSplit[1];
        String chainId = filenameSplit[2].split("-")[0];
        try {
            String matchEcNumber = SiftsMappingAnnotator.getEnzymeLineForPdbId(pdbId, chainId)
                    .map(split -> split[3])
                    .orElse("X.X.X.X");
            String classLabel = inferClassLabel(originalEcNumber, matchEcNumber);

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();
            int[] residueNumbers = IntStream.range(2, filenameSplit.length)
                    .mapToObj(index -> filenameSplit[index])
                    .map(string -> string.split("-")[1])
                    .mapToInt(Integer::valueOf)
                    .toArray();
            System.out.println(Arrays.toString(residueNumbers));
            List<Group> aminoAcids = chain.select()
                    .residueNumber(residueNumbers)
                    .asFilteredGroups()
                    .collect(Collectors.toList());
            int size = aminoAcids.size();

            ACCESSIBLE_SURFACE_AREA_CALCULATOR.process(structure);
            DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE.process(structure);
            ENERGY_PROFILE_CALCULATOR.process(structure);
            LOOP_FRACTION_CALCULATOR.process(structure);
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(structure);

            double[] rasa = aminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(AccessibleSurfaceArea.class))
                    .mapToDouble(AccessibleSurfaceArea::getRelativeAccessibleSurfaceArea)
                    .toArray();
            double[] energy = aminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(EnergyProfile.class))
                    .mapToDouble(EnergyProfile::getSolvationEnergy)
                    .toArray();
            double[] loop = aminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(AccessibleSurfaceArea.class))
                    .mapToDouble(AccessibleSurfaceArea::getRelativeAccessibleSurfaceArea)
                    .toArray();
            int[] hydrogenBonds = aminoAcids.stream()
                    .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                    .mapToInt(container -> container.getHydrogenBonds().size())
                    .toArray();

            return Optional.of(motifId + "," +
                    IntStream.range(0, size).mapToObj(index -> StandardFormat.format(rasa[index])).collect(Collectors.joining(",")) + "," +
                    IntStream.range(0, size).mapToObj(index -> StandardFormat.format(energy[index])).collect(Collectors.joining(",")) + "," +
                    IntStream.range(0, size).mapToObj(index -> StandardFormat.format(loop[index])).collect(Collectors.joining(",")) + "," +
                    IntStream.range(0, size).mapToObj(index -> StandardFormat.format(hydrogenBonds[index])).collect(Collectors.joining(",")) + "," +
                    classLabel);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static String inferClassLabel(String originalEcNumber, String matchEcNumber) {
        String[] originalEcNumberSplit = originalEcNumber.split("\\.");
        String[] matchEcNumberSplit = matchEcNumber.split("\\.");
        if(!originalEcNumberSplit[0].equals(matchEcNumberSplit[0])) {
            return "non-functional";
        }
        if(!originalEcNumberSplit[1].equals(matchEcNumberSplit[1])) {
            return "non-functional";
        }
        if(!originalEcNumberSplit[2].equals(matchEcNumberSplit[2])) {
            return "non-functional";
        }
        return "functional";
    }
}

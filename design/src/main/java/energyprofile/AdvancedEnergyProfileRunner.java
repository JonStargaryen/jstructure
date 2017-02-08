package energyprofile;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileAligner;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.parser.CIFParser;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;

/**
 * Provides more command line options to the user.
 * Created by bittrich on 2/6/17.
 */
public class AdvancedEnergyProfileRunner {
    public static void main(String... args) {
        if(args.length == 0 || args.length > 3) {
            System.out.println(printUsage());
            throw new IllegalArgumentException();
        }

        builder()
                .first(args[0])
                .output(args[args.length - 1])
                .second(args[1])
                .execute();
    }

    private static String printUsage() {
        return "too few command-line arguments" + System.lineSeparator() +
                "the jar provides 3 modes to handle energy profiles, the last argument is always the output file's destination:" + System.lineSeparator() +
                "\tcompute an energy profile for a given structure" + System.lineSeparator() +
                "\tpredict an energy profile for a given sequence or structure (solely using its sequence)" + System.lineSeparator() +
                "\talign two energy profiles" + System.lineSeparator() +
                "/path/to/file.pdb /path/to/output.ep - compute the energy profile and write the output to the specified path" + System.lineSeparator() +
                "MILIAQERTW... /path/to/output.ep - predict the energy profile for the given sequence" + System.lineSeparator() +
                "/path/to/file1.pdb /path/to/file2.pdb /path/to/output.epa - compute and align both energy profiles, write results" + System.lineSeparator() +
                "/path/to/file.pdb MILIAQERTW... /path/to/output.epa - compute/predict energy profile and align";
    }

    private static Builder builder() {
        return new Builder();
    }

    static class Builder {
        private static final AbstractFeatureProvider energyProfileCalculator = FeatureProviderRegistry.resolveAnnotator(EnergyProfileCalculator.SOLVATION_ENERGY);
        private static final AbstractFeatureProvider energyProfilePredictor = FeatureProviderRegistry.resolvePredictor(EnergyProfileCalculator.SOLVATION_ENERGY);
        private static final EnergyProfileAligner energyProfileAligner = new EnergyProfileAligner();
        private ProteinContainer proteinContainer1;
        private ProteinContainer proteinContainer2;
        private Path outputPath;
        private boolean alignmentMode;

        Builder output(String argument) {
            outputPath = Paths.get(argument);
            return this;
        }

        Builder first(String argument) {
            proteinContainer1 = handleProteinReference(argument);
            return this;
        }

        Builder second(String argument) {
            // if this argument already describes the output path, then do nothing
            if(Paths.get(argument).equals(outputPath)) {
                return this;
            }

            proteinContainer2 = handleProteinReference(argument);
            alignmentMode = true;
            return this;
        }

        void execute() {
            String output;
            if(alignmentMode) {
                double distanceScore = energyProfileAligner.align(proteinContainer1.getEnergyProfile(), proteinContainer2.getEnergyProfile());
                output = "energy profile alignment between " + proteinContainer1.protein.getIdentifier() + " and " + proteinContainer2.protein.getIdentifier() + System.lineSeparator() +
                        "distance score: " + distanceScore + System.lineSeparator() +
                        "continuous values between 0 and 5, 0 = optimal alignment, 5 = worst possible alignment score";
            } else {
                output = composeOutput(proteinContainer1.protein);
            }

            try {
                Files.write(outputPath, output.getBytes());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        private ProteinContainer handleProteinReference(String argument) {
            Path potentialPath = Paths.get(argument);
            // argument refers to a file
            if(Files.exists(potentialPath)) {
                return new ProteinContainer(potentialPath);
            }

            // argument refers to a protein sequence
            if(argument.matches("^[ACDEFGHIKLMNPQRSTVWY]+$")) {
                return new ProteinContainer(argument);
            }

            throw new IllegalArgumentException("provided argument '" + argument + "' does neither match a file nor is a protein sequence");
        }
    }

    static class ProteinContainer {
        Path pathToPdbFile;
        Protein protein;

        List<Double> getEnergyProfile() {
            return protein.aminoAcids()
                    .map(this::getSolvationEnergy)
                    .collect(Collectors.toList());
        }

        private double getSolvationEnergy(Group group) {
            return group.getFeatureAsDouble(EnergyProfileCalculator.SOLVATION_ENERGY);
        }

        ProteinContainer(Path pathToPdbFile) {
            this.pathToPdbFile = pathToPdbFile;
            this.protein = ProteinParser.parsePDBFile(pathToPdbFile);
            Builder.energyProfileCalculator.process(protein);
        }

        ProteinContainer(String sequence) {
            this.protein = new Protein();
            this.protein.setName("UNKNOWN PROTEIN");
            Chain chain = new Chain("A");
            this.protein.addChain(chain);
            for(int resNum = 1; resNum <= sequence.length(); resNum++) {
                String aminoAcidName = String.valueOf(sequence.charAt(resNum - 1));
                AminoAcidFamily aminoAcidFamily = AminoAcidFamily.valueOfIgnoreCase(aminoAcidName).orElseThrow(() -> new IllegalArgumentException("unknown amino acid '" + aminoAcidName + "'"));
                Group group = new Group(aminoAcidFamily.getThreeLetterCode(), resNum, CIFParser.parseLigandInformation(aminoAcidFamily.getThreeLetterCode()), false);
                chain.addGroup(group);
            }
            Builder.energyProfilePredictor.process(protein);
        }
    }

    private static final String DELIMITER = "\t";
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.US));

    static String composeOutput(Protein protein) {
        return protein.aminoAcids()
                .map(AdvancedEnergyProfileRunner::composeOutputLine)
                .collect(Collectors.joining(System.lineSeparator(),
                        "ID" + DELIMITER + protein.getIdentifier() + System.lineSeparator() +
                                "chain" + DELIMITER + "resNum" + DELIMITER + "resName" + DELIMITER + "energy" + System.lineSeparator(),
                        ""));
    }

    private static String composeOutputLine(Group group) {
        return group.getParentChain().getIdentifier() + DELIMITER + group.getResidueNumber() + DELIMITER +
                group.getThreeLetterCode() + DELIMITER + DECIMAL_FORMAT.format(group.getFeatureAsDouble(EnergyProfileCalculator.SOLVATION_ENERGY));
    }
}

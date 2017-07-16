package studies.poznan;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileAligner;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ProteinIdentifier;

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
                "\tcalculate an energy profile for a given structure" + System.lineSeparator() +
                "\tpredict an energy profile for a given sequence or structure (solely using its sequence)" + System.lineSeparator() +
                "\talign two energy profiles" + System.lineSeparator() +
                "/path/to/file.pdb /path/to/output.ep - calculate the energy profile and write the output to the specified path" + System.lineSeparator() +
                "MILIAQERTW... /path/to/output.ep - predict the energy profile for the given sequence" + System.lineSeparator() +
                "/path/to/file1.pdb /path/to/file2.pdb /path/to/output.epa - calculate and align both energy profiles, write results" + System.lineSeparator() +
                "/path/to/file.pdb MILIAQERTW... /path/to/output.epa - calculate/predict energy profile and align";
    }

    private static Builder builder() {
        return new Builder();
    }

    static class Builder {
        private static final AbstractFeatureProvider energyProfileCalculator = FeatureProviderRegistry.resolveAnnotator(EnergyProfile.class);
        private static final AbstractFeatureProvider energyProfilePredictor = FeatureProviderRegistry.resolvePredictor(EnergyProfile.class);
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
                double distanceScore = energyProfileAligner.align(proteinContainer1.getEnergyProfile(),
                        proteinContainer2.getEnergyProfile());
                output = "energy profile alignment between " + proteinContainer1.protein.getIdentifier() + " and " +
                        proteinContainer2.protein.getIdentifier() + System.lineSeparator() +
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
        Structure protein;

        List<Double> getEnergyProfile() {
            return protein.aminoAcids()
                    .map(this::getSolvationEnergy)
                    .collect(Collectors.toList());
        }

        private double getSolvationEnergy(Group group) {
            return group.getFeatureContainer().getFeature(EnergyProfile.class).getSolvationEnergy();
        }

        ProteinContainer(Path pathToPdbFile) {
            this.pathToPdbFile = pathToPdbFile;
            this.protein = StructureParser.source(pathToPdbFile).parse();
            Builder.energyProfileCalculator.process(protein);
        }

        ProteinContainer(String sequence) {
            this.protein = new Structure(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER);
            Chain chain = new Chain(IdentifierFactory.createChainIdentifier(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER, "A"));
            this.protein.addChain(chain);
            for(int resNum = 1; resNum <= sequence.length(); resNum++) {
                String aminoAcidName = String.valueOf(sequence.charAt(resNum - 1));
                AminoAcid.Family aminoAcid = AminoAcid.Family.resolveOneLetterCode(aminoAcidName);
                Group group = new Group(aminoAcid.getGroupPrototype(),
                        IdentifierFactory.createResidueIdentifier(resNum),
                        false);
                chain.addGroup(group);
            }
            Builder.energyProfilePredictor.process(protein);
        }
    }

    private static final String DELIMITER = "\t";
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000",
            DecimalFormatSymbols.getInstance(Locale.US));

    static String composeOutput(Structure protein) {
        return protein.aminoAcids()
                .map(AdvancedEnergyProfileRunner::composeOutputLine)
                .collect(Collectors.joining(System.lineSeparator(),
                        "ID" + DELIMITER + protein.getIdentifier() + System.lineSeparator() +
                                "chain" + DELIMITER + "resNum" + DELIMITER + "resName" + DELIMITER + "energy" +
                                System.lineSeparator(),
                        ""));
    }

    private static String composeOutputLine(Group group) {
        return group.getParentChain().getIdentifier() + DELIMITER + group.getResidueIdentifier() + DELIMITER +
                group.getThreeLetterCode() + DELIMITER +
                DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(EnergyProfile.class).getSolvationEnergy());
    }
}

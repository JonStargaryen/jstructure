package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.stream.Collectors;

/**
 * A wrapper for the energy profile calculation so it can be accessed as CLI.
 * Created by S on 16.01.2017.
 */
public class EnergyProfileCalculatorRunner {
    private static final String DELIMITER = "\t";
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000",DecimalFormatSymbols.getInstance(Locale.US));

    public static void main(String... args) {
        if(args.length != 1 && args.length != 2) {
            // print usage
            System.err.println("java -jar energy.jar /path/to/file.pdb [/path/to/output/file/if/desired.ep2]");
            throw new IllegalArgumentException("too " + (args.length < 1 ? "few" : "many") + " arguments provided");
        }

        String inputPath = args[0];
        boolean outputRequested = false;
        String outputPath = null;

        if(args.length == 2) {
            outputRequested = true;
            outputPath = args[1];
        }

        try {
            // read and parse protein
            Protein protein = ProteinParser.parsePDBFile(inputPath);

            // resolve feature provider by name
            AbstractFeatureProvider featureProvider = FeatureProviderRegistry.getInstance().resolve(EnergyProfileCalculator.SOLVATION_ENERGY);

            // compute energy profile
            featureProvider.process(protein);

            // compose output
            String output = protein.aminoAcids()
                    .map(EnergyProfileCalculatorRunner::composeOutputLine)
                    .collect(Collectors.joining(System.lineSeparator(),
                            "ID" + DELIMITER + protein.getIdentifier() + System.lineSeparator() +
                            "chain" + DELIMITER + "resNum" + DELIMITER + "resName" + DELIMITER + "energy" + System.lineSeparator(),
                            ""));

            System.out.println(output);

            if(outputRequested) {
                Files.write(Paths.get(outputPath), output.getBytes());
            }
        } catch (UncheckedIOException e) {
            System.err.println("failed to parse PDB file located at " + inputPath);
            throw e;
        } catch (IOException e) {
            System.err.println("writing the output file failed - destination " + outputPath);
            throw new UncheckedIOException(e);
        }
    }

    private static String composeOutputLine(Group group) {
        return group.getParentChain().getIdentifier() + DELIMITER + group.getResidueNumber() + DELIMITER +
                group.getThreeLetterCode() + DELIMITER + DECIMAL_FORMAT.format(group.getFeatureAsDouble(EnergyProfileCalculator.SOLVATION_ENERGY));
    }
}

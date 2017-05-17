package studies.poznan;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.apache.log4j.BasicConfigurator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * A wrapper for the energy profile calculation so it can be accessed as CLI.
 * Created by S on 16.01.2017.
 */
public class SimpleEnergyProfileRunner {
    private static final Logger logger = LoggerFactory.getLogger(SimpleEnergyProfileRunner.class);
    private static final AbstractFeatureProvider featureProvider = FeatureProviderRegistry.resolve(EnergyProfile.class);

    static {
        BasicConfigurator.configure();
    }

    public static void main(String... args) {
        if(args.length != 1 && args.length != 2) {
            // print usage
            System.err.println("usage: java -jar energy.jar /path/to/file.pdb [/path/to/output/file/if/desired.ep]");
            throw new IllegalArgumentException("too " + (args.length < 1 ? "few" : "many") + " arguments provided");
        }

        String inputPath = args[0];
        boolean outputRequested = false;
        String outputPath = null;

        if(args.length == 2) {
            outputRequested = true;
            outputPath = args[1];
        }

        logger.info("computing energy profile for {}", inputPath);
        logger.info("writing output to {}", outputRequested ? outputPath : "console");

        try {
            // read and parse protein
            Protein protein = ProteinParser.source(inputPath).parse();

            // compute energy profile
            featureProvider.process(protein);

            // compose output
            String output = AdvancedEnergyProfileRunner.composeOutput(protein);

            if(outputRequested) {
                Files.write(Paths.get(outputPath), output.getBytes());
            } else {
                System.out.println(output);
            }
        } catch (UncheckedIOException e) {
            System.err.println("failed to parse PDB file located at " + inputPath);
            throw e;
        } catch (IOException e) {
            System.err.println("writing the output file failed - destination " + outputPath);
            throw new UncheckedIOException(e);
        }
    }
}

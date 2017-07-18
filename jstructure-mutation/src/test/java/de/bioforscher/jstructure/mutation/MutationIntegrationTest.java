package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.impl.MutationEffectPredictionServiceImpl;
import de.bioforscher.jstructure.mutation.impl.MutationFeatureVectorImpl;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Test the complex interactions of this module's classes. Should run the test datasets.
 * Created by bittrich on 7/17/17.
 */
public class MutationIntegrationTest {
    private MutationEffectPredictionService mutationEffectPredictionService;

    @Before
    public void setup() {
        mutationEffectPredictionService = new MutationEffectPredictionServiceImpl();
    }

    @Test
    public void shouldComposeArffForSchaefer2012Stability() throws IOException {
        // group entries by their pdbId so everything has to be computed only once for changing positions or chains
        Map<String, List<String[]>> sequenceMap = TestUtils.getResourceAsStream("schaefer2012/mutants_stability.txt")
                .filter(line -> !line.startsWith("ID,"))
                .map(line -> line.split(","))
                .collect(Collectors.groupingBy(split -> split[0]));

        // write header if needed, otherwise append file
        Path outputPath = Paths.get("/home/bittrich/schaefer2012-structure.arff");
        boolean writeHeader = !Files.exists(outputPath);
        FileWriter fileWriter = new FileWriter(outputPath.toFile(), true);
        if (writeHeader) {
            String header = getHeaderSchaefer2012Stability();
            fileWriter.append(header);
            System.out.println("writing header:" + System.lineSeparator() + System.lineSeparator() + header);
        }

        // handle individual sequence bins
        sequenceMap.entrySet().forEach(group -> handleMutantLineGroup(group, fileWriter));
    }

    private void handleMutantLineGroup(Map.Entry<String, List<String[]>> entry, FileWriter fileWriter) {
        String jobName = entry.getKey();
        String pdbId = jobName.substring(0, 4);
        String chainId = jobName.substring(4);

        // dropping other chains for consistent computation
        Structure protein = StructureParser.source(pdbId)
                .minimalParsing(true)
                .parse()
                .select()
                .chainId(chainId)
                .asIsolatedStructure();
        Chain chain = protein.select()
                .chainName(chainId)
                .asChain();

        // create global mutation job instance shared by all entries of this bin
        MutationJob mutationJob = mutationEffectPredictionService.createMutationJob(jobName, chain);

        String partialOutput = entry.getValue().stream()
                .map(mutantLine -> handleMutantLine(mutationJob, mutantLine))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .peek(System.out::println)
                .collect(Collectors.joining(System.lineSeparator()));

        try {
            fileWriter.write(partialOutput);
            fileWriter.flush();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }


    private Optional<String> handleMutantLine(MutationJob mutationJob, String[] mutantLine) {
        try {
            ResidueIdentifier residueIdentifier = IdentifierFactory.createResidueIdentifier(Integer.valueOf(mutantLine[1]));
            AminoAcid.Family targetAminoAcid = AminoAcid.Family.resolveOneLetterCode(mutantLine[2]);
            double ddG = Double.valueOf(mutantLine[3]);
            String classLabel = Math.abs(ddG) > 1.0 ? "effect" : "tolerated";

            // create feature vector for this particularized mutation
            MutationFeatureVector mutationEffectVector = mutationEffectPredictionService.createMutationFeatureVector(mutationJob,
                    mutationJob.getReferenceChain().getChainIdentifier(),
                    residueIdentifier,
                    targetAminoAcid);

            return Optional.of(mutationEffectVector.toString() + "," + classLabel);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    public String getHeaderSchaefer2012Stability() {
        return "@RELATION mutation" + System.lineSeparator() +
                "@ATTRIBUTE chainId        STRING" + System.lineSeparator() +
                "@ATTRIBUTE mutation       STRING" + System.lineSeparator() +
                MutationFeatureVectorImpl.toPartialHeader() + System.lineSeparator() +
                "@ATTRIBUTE class          {effect,tolerated}" + System.lineSeparator() +
                "@DATA" + System.lineSeparator();
    }
}
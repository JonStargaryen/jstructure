package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.align.impl.TMAlignService;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.si.ConfoldServiceWorker;
import de.bioforscher.jstructure.si.explorer.DataSource;
import de.bioforscher.jstructure.si.explorer.ExplorerChain;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class A03C_MissingStrategyRunner {
    private static final Logger logger = LoggerFactory.getLogger(A03C_MissingStrategyRunner.class);
    private static final double DEFAULT_COVERAGE = 0.3;
    private static final int REDUNDANCY = 10;
    private static final Path OUTPUT_PATH = Paths.get("/home/sb/reconstruction.csv");
    private static final TMAlignService TM_ALIGN_SERVICE = TMAlignService.getInstance();

    private static FileWriter fileWriter;
    private static ExecutorService executorService;

    private static final List<String> MISSING_COMBINATIONS = Pattern.compile("\n").splitAsStream("STF0001,IgnoreWorst\n" +
            "STF0003,IgnoreWorst\n" +
            "STF0004,IgnoreWorst\n" +
            "STF0006,IgnoreWorst\n" +
            "STF0008,IgnoreWorst\n" +
            "STF0009,IgnoreWorst\n" +
            "STF0010,IgnoreWorst\n" +
            "STF0011,IgnoreWorst\n" +
            "STF0012,IgnoreWorst\n" +
            "STF0013,IgnoreWorst\n" +
            "STF0014,IgnoreWorst\n" +
            "STF0015,IgnoreWorst\n" +
            "STF0016,IgnoreWorst\n" +
            "STF0018,IgnoreWorst\n" +
            "STF0019,BestByAverage\n" +
            "STF0019,IgnoreBest\n" +
            "STF0019,Random\n" +
            "STF0019,WorstByAverage\n" +
            "STF0019,NonNative\n" +
            "STF0019,All\n" +
            "STF0019,IgnoreWorst\n" +
            "STF0020,BestByAverage\n" +
            "STF0020,IgnoreBest\n" +
            "STF0020,Random\n" +
            "STF0020,WorstByAverage\n" +
            "STF0020,NonNative\n" +
            "STF0020,All\n" +
            "STF0020,IgnoreWorst\n" +
            "STF0021,IgnoreWorst\n" +
            "STF0023,IgnoreWorst\n" +
            "STF0024,IgnoreWorst\n" +
            "STF0025,IgnoreWorst\n" +
            "STF0026,IgnoreWorst\n" +
            "STF0028,IgnoreWorst\n" +
            "STF0037,IgnoreWorst\n" +
            "STF0038,IgnoreWorst\n" +
            "STF0040,IgnoreWorst\n" +
            "STF0042,IgnoreWorst\n" +
            "STF0043,IgnoreWorst\n" +
            "STF0044,IgnoreWorst\n" +
            "STF0045,IgnoreWorst\n" +
            "STF0046,BestByAverage\n" +
            "STF0046,IgnoreBest\n" +
            "STF0046,Random\n" +
            "STF0046,WorstByAverage\n" +
            "STF0046,NonNative\n" +
            "STF0046,All\n" +
            "STF0046,IgnoreWorst")
            .collect(Collectors.toList());

    public static void main(String[] args) throws IOException {
            fileWriter = new FileWriter(OUTPUT_PATH.toFile());
            fileWriter.write("id,strategy,rmsd" + System.lineSeparator());
            executorService = Executors.newFixedThreadPool(16);

        DataSource.getInstance()
                .start2FoldChains()
                .forEach(A03C_MissingStrategyRunner::handleChain);
    }

    private static void handleChain(ExplorerChain explorerChain) {
        logger.info("handling chain {}",
                explorerChain.getStfId());

        try {
            Chain nativeChain = explorerChain.getChain();
            Path nativeChainPath = Files.createTempFile("nativechain-", ".pdb");
            Files.write(nativeChainPath, nativeChain.getPdbRepresentation().getBytes());

            List<ContactStructuralInformation> contactStructuralInformation = explorerChain.getContacts();
            int numberOfNativeContacts = contactStructuralInformation.size();
            int numberOfContactsToSelect = (int) (numberOfNativeContacts * DEFAULT_COVERAGE);

            List<ReconstructionContactMap> contactMaps = Stream.of(A03_ReconstructByVariousStrategy.ReconstructionStrategyDefinition.values())
                    .map(A03_ReconstructByVariousStrategy.ReconstructionStrategyDefinition::getReconstructionStrategy)
                    // filter for missing combinations
                    .filter(strategy -> MISSING_COMBINATIONS.stream()
                            .filter(combination -> combination.startsWith(explorerChain.getStfId()))
                            .anyMatch(combination -> combination.endsWith(strategy.getClass().getSimpleName())))
                    .flatMap(reconstructionStrategy -> IntStream.range(0, REDUNDANCY)
                            .mapToObj(i -> {
                                ReconstructionContactMap contactMap = reconstructionStrategy.composeReconstructionContactMap(nativeChain,
                                        contactStructuralInformation,
                                        numberOfContactsToSelect);

                                contactMap.setName(reconstructionStrategy.getClass().getSimpleName() + "-" + (i + 1));

                                return contactMap;
                            }))
                    .collect(Collectors.toList());

            // nothing to do
            if(contactMaps.isEmpty()) {
                return;
            }

            Map<String, List<Future<List<Chain>>>> reconstructionFutures = new HashMap<>();
            for (ReconstructionContactMap contactMap : contactMaps) {
                String name = contactMap.getName().split("-")[0];
                logger.info("handling contact map definition {}",
                        name);

                if(!reconstructionFutures.containsKey(name)) {
                    reconstructionFutures.put(name, new ArrayList<>());
                }

                List<Future<List<Chain>>> bin = reconstructionFutures.get(name);

                bin.add(executorService.submit(new ConfoldServiceWorker("/home/sb/programs/confold_v1.0/confold.pl",
                        contactMap.getSequence(),
                        contactMap.getSecondaryStructureElements(),
                        contactMap.getCaspRRRepresentation())));
            }

            for (Map.Entry<String, List<Future<List<Chain>>>> reconstructionFuture : reconstructionFutures.entrySet()) {
                try {
                    String name = reconstructionFuture.getKey();
                    List<Chain> reconstructions = reconstructionFuture.getValue()
                            .stream()
                            .map(future -> {
                                try {
                                    return future.get();
                                } catch (Exception e) {
                                    throw new ComputationException(e);
                                }
                            })
                            .flatMap(Collection::stream)
                            .collect(Collectors.toList());
                    List<TMAlignAlignmentResult> alignmentResults = new ArrayList<>();
                    List<Path> tmpFiles = new ArrayList<>();

                    if (reconstructions.isEmpty()) {
                        throw new ComputationException("reconstruction did not yield any reconstructs");
                    }

                    for (Chain reconstructedChain : reconstructions) {
                        Path reconstructPath = Files.createTempFile("confoldservice-recon", ".pdb");
                        tmpFiles.add(reconstructPath);
                        Files.write(reconstructPath, reconstructedChain.getPdbRepresentation().getBytes());
                        alignmentResults.add(TM_ALIGN_SERVICE.process(new String[]{
                                "/home/sb/programs/tmalign",
                                nativeChainPath.toFile().getAbsolutePath(),
                                reconstructPath.toFile().getAbsolutePath()
                        }));
                    }

                    if (alignmentResults.isEmpty()) {
                        throw new ComputationException("tmalign did not yield any alignments");
                    }

                    for(TMAlignAlignmentResult alignmentResult : alignmentResults) {
                        double rmsd = alignmentResult.getRootMeanSquareDeviation().getScore();
                        String line = explorerChain.getStfId() + "," + name + "," + rmsd;
                        logger.info(line);
                        fileWriter.write(line + System.lineSeparator());
                        fileWriter.flush();
                    }

                    // cleanup
                    for (Path tmpFile : tmpFiles) {
                        Files.delete(tmpFile);
                    }
                } catch (IOException e) {
                    throw new ComputationException(e);
                }
            }
        } catch (IOException e) {
            throw new ComputationException(e);
        }
    }
}

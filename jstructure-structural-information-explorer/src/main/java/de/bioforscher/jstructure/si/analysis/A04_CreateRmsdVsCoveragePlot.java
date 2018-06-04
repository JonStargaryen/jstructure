package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.align.impl.TMAlignService;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.ConfoldServiceWorker;
import de.bioforscher.jstructure.si.explorer.DataSource;
import de.bioforscher.jstructure.si.explorer.ExplorerChain;
import de.bioforscher.jstructure.si.model.ReconstructionResult;
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
import java.util.stream.Collectors;

public class A04_CreateRmsdVsCoveragePlot {
    private static final Logger logger = LoggerFactory.getLogger(A04_CreateRmsdVsCoveragePlot.class);
    private static final int REDUNDANCY = 10;
    private static final Path OUTPUT_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/si/reconstruction-perc.csv");
    private static final TMAlignService TM_ALIGN_SERVICE = TMAlignService.getInstance();

    private static FileWriter fileWriter;
    private static ExecutorService executorService;

    public static void main(String[] args) throws IOException {
        boolean fileExisted = Files.exists(OUTPUT_PATH);
        fileWriter = new FileWriter(OUTPUT_PATH.toFile(), true);
        if(!fileExisted) {
            fileWriter.write("id,coverage,rmsd" + System.lineSeparator());
        }
        executorService = Executors.newFixedThreadPool(8);

        DataSource.getInstance()
                .start2FoldChains()
                // skip failing computation STF0006
                .filter(explorerChain -> !explorerChain.getStfId().equals("STF0006") && !explorerChain.getStfId().equals("STF0044"))
                .forEach(A04_CreateRmsdVsCoveragePlot::handleChain);
    }

    private static void handleChain(ExplorerChain explorerChain) {
        logger.info("handling chain {}",
                explorerChain.getStfId());

        try {
            Chain nativeChain = explorerChain.getChain();
            Path nativeChainPath = Files.createTempFile("nativechain-", ".pdb");
            Files.write(nativeChainPath, nativeChain.getPdbRepresentation().getBytes());

            ReconstructionContactMap nativeContactMap = ReconstructionContactMap.createReconstructionContactMap(nativeChain,
                    ContactDefinitionFactory.createAlphaCarbonContactDefinition(8.0));
            List<AminoAcid> aminoAcids = nativeChain.getAminoAcids();
            List<Pair<AminoAcid, AminoAcid>> contacts = nativeContactMap.getLongRangeContacts();
            int numberNativeLongRangeContacts = contacts.size();
            List<ReconstructionContactMap> reconstructionContactMaps = new ArrayList<>();

            for(int coverage = 5; coverage <= 100; coverage = coverage + 5) {
                int numberOfContactsToSelect = (int) Math.round(0.01 * coverage * numberNativeLongRangeContacts);
                for(int run = 0; run < REDUNDANCY; run++) {
                    Collections.shuffle(contacts);
                    List<Pair<AminoAcid, AminoAcid>> selectedContacts = contacts.subList(0, numberOfContactsToSelect);
                    ReconstructionContactMap contactMap = new ReconstructionContactMap(aminoAcids, selectedContacts, nativeContactMap.getContactDefinition());
                    contactMap.setName("p" + coverage + "-" + (run + 1));
                    reconstructionContactMaps.add(contactMap);
                }
            }

            Map<String, List<Future<ReconstructionResult>>> reconstructionFutures = new HashMap<>();
            for (ReconstructionContactMap contactMap : reconstructionContactMaps) {
                String name = contactMap.getName().split("-")[0];
                logger.info("handling contact map with coverage {}",
                        name);

                if(!reconstructionFutures.containsKey(name)) {
                    reconstructionFutures.put(name, new ArrayList<>());
                }

                if(Files.lines(OUTPUT_PATH)
                        .anyMatch(line -> line.startsWith(explorerChain.getStfId()) && line.split(",")[1].equals(name))) {
                    continue;
                }

                List<Future<ReconstructionResult>> bin = reconstructionFutures.get(name);

                bin.add(executorService.submit(new ConfoldServiceWorker("/home/bittrich/programs/confold_v1.0/confold.pl",
                        contactMap.getSequence(),
                        contactMap.getSecondaryStructureElements(),
                        contactMap.getCaspRRRepresentation(),
                        nativeContactMap.getConfoldRRType())));
            }

            for (Map.Entry<String, List<Future<ReconstructionResult>>> reconstructionFuture : reconstructionFutures.entrySet()) {
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
                            .map(ReconstructionResult::getChains)
                            .flatMap(Collection::stream)
                            .collect(Collectors.toList());
                    List<TMAlignAlignmentResult> alignmentResults = new ArrayList<>();
                    List<Path> tmpFiles = new ArrayList<>();

                    if (reconstructions.isEmpty()) {
                        // bin already processed - skipping this may also result in skipping upon failed computation
                        continue;
                    }

                    for (Chain reconstructedChain : reconstructions) {
                        Path reconstructPath = Files.createTempFile("confoldservice-recon", ".pdb");
                        tmpFiles.add(reconstructPath);
                        Files.write(reconstructPath, reconstructedChain.getPdbRepresentation().getBytes());
                        alignmentResults.add(TM_ALIGN_SERVICE.process(new String[] {
                                "/home/bittrich/programs/tmalign/tmalign",
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
        } catch (IOException | AlignmentException e) {
            throw new ComputationException(e);
        }
    }
}

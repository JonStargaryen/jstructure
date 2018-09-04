package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.ConfoldServiceWorker;
import de.bioforscher.jstructure.si.explorer.ExplorerChain;
import de.bioforscher.jstructure.si.model.ReconstructionResult;
import de.bioforscher.testutil.TestUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class A04A_CreatePyMolRendering {
    private static final Logger logger = LoggerFactory.getLogger(A04A_CreatePyMolRendering.class);
    private static final int REDUNDANCY = 10;
    private static final Path OUTPUT_PATH = Paths.get("/home/sb/reconstructions/");
    private static ExecutorService executorService;

    public static void main(String[] args) throws IOException {
        boolean fileExisted = Files.exists(OUTPUT_PATH);
        executorService = Executors.newFixedThreadPool(16);

        TestUtils.getResourceAsStream("data/efr.list")
                .map(line -> line.split(";"))
                .filter(split -> split[0].equals("STF0023"))
                .map(ExplorerChain::new)
                .forEach(A04A_CreatePyMolRendering::handleChain);
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

            IntStream.of(5, 30, 100).forEach(coverage -> {
                int numberOfContactsToSelect = (int) Math.round(0.01 * coverage * numberNativeLongRangeContacts);
                for(int run = 0; run < REDUNDANCY; run++) {
                    Collections.shuffle(contacts);
                    List<Pair<AminoAcid, AminoAcid>> selectedContacts = contacts.subList(0, numberOfContactsToSelect);
                    ReconstructionContactMap contactMap = new ReconstructionContactMap(aminoAcids, selectedContacts, nativeContactMap.getContactDefinition());
                    contactMap.setName("p" + coverage + "-" + (run + 1));
                    reconstructionContactMaps.add(contactMap);
                }
            });

            Map<String, List<Future<ReconstructionResult>>> reconstructionFutures = new HashMap<>();
            for (ReconstructionContactMap contactMap : reconstructionContactMaps) {
                String name = contactMap.getName().split("-")[0];
                logger.info("handling contact map with coverage {}",
                        name);

                if(!reconstructionFutures.containsKey(name)) {
                    reconstructionFutures.put(name, new ArrayList<>());
                }

                List<Future<ReconstructionResult>> bin = reconstructionFutures.get(name);

                bin.add(executorService.submit(new ConfoldServiceWorker("/home/sb/programs/confold_v1.0/confold.pl",
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

                    for (Chain reconstructedChain : reconstructions) {
                        Files.write(OUTPUT_PATH.resolve(name + ".pdb"), reconstructedChain.getPdbRepresentation().getBytes());
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

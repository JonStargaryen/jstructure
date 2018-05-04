package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.align.impl.TMAlignService;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
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
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class A03_ReconstructByVariousStrategy {
    private static final Logger logger = LoggerFactory.getLogger(A03_ReconstructByVariousStrategy.class);
    private static final List<String> INTERESTING_STRUCTURES = Stream.of("1lqv",
            "1mjc",
            "1ten",
            "2acy",
            "3chy",
            "1rc2",
            "STF0021",
            "STF0023",
            "STF0042")
            .collect(Collectors.toList());
    private static final double DEFAULT_COVERAGE = 0.3;
    private static final int REDUNDANCY = 10;
    private static final Path OUTPUT_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/si/reconstruction.csv");
    private static final TMAlignService TM_ALIGN_SERVICE = TMAlignService.getInstance();

    private static FileWriter fileWriter;
    private static ExecutorService executorService;

    public static void main(String[] args) throws IOException {
        fileWriter = new FileWriter(OUTPUT_PATH.toFile());
        fileWriter.write("id,strategy,rmsd" + System.lineSeparator());
        executorService = Executors.newFixedThreadPool(8);

        DataSource.getInstance()
                .chains()
                .filter(explorerChain -> INTERESTING_STRUCTURES.stream()
                        .noneMatch(id -> explorerChain.getStfId().startsWith(id)))
                .forEach(A03_ReconstructByVariousStrategy::handleChain);
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

            List<ReconstructionContactMap> contactMaps = Stream.of(ReconstructionStrategyDefinition.values())
                    .map(ReconstructionStrategyDefinition::getReconstructionStrategy)
                    .flatMap(reconstructionStrategy -> IntStream.range(0, REDUNDANCY)
                            .mapToObj(i -> {
                                ReconstructionContactMap contactMap = reconstructionStrategy.composeReconstructionContactMap(nativeChain,
                                        contactStructuralInformation,
                                        numberOfContactsToSelect);

                                contactMap.setName(reconstructionStrategy.getClass().getSimpleName() + "-" + (i + 1));

                                return contactMap;
                            }))
                    .collect(Collectors.toList());

            Map<String, List<Future<List<Chain>>>> reconstructionFutures = new HashMap<>();
            for (ReconstructionContactMap contactMap : contactMaps) {
                String name = contactMap.getName().split("-")[0];
                logger.info("handling contact map definition {}",
                        name);

                if(!reconstructionFutures.containsKey(name)) {
                    reconstructionFutures.put(name, new ArrayList<>());
                }

                List<Future<List<Chain>>> bin = reconstructionFutures.get(name);

                bin.add(executorService.submit(new ConfoldServiceWorker("/home/bittrich/programs/confold_v1.0/confold.pl",
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
        } catch (IOException e) {
            throw new ComputationException(e);
        }
    }

    interface ReconstructionStrategy {
        List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                     List<ContactStructuralInformation> contactStructuralInformation,
                                                     int numberOfContacts);

        default ReconstructionContactMap composeReconstructionContactMap(Chain chain,
                                                                         List<ContactStructuralInformation> contactStructuralInformation,
                                                                         int numberOfContacts) {
            return new ReconstructionContactMap(chain.getAminoAcids(),
                    selectContacts(chain, contactStructuralInformation, numberOfContacts)
                            .stream()
                            .map(contact -> new Pair<>(chain.select()
                                            .residueNumber(contact.getLeft())
                                            .asAminoAcid(),
                                    chain.select()
                                            .residueNumber(contact.getRight())
                                            .asAminoAcid()))
                            .collect(Collectors.toList()));
        }

        String getName();
    }

    enum ReconstructionStrategyDefinition {
        BEST_BY_AVERAGE(new BestByAverage()),
        IGNORE_BEST(new IgnoreBest()),
        RANDOM(new Random()),
        WORST_BY_AVERAGE(new WorstByAverage()),
        NON_NATIVE(new NonNative());

        private ReconstructionStrategy reconstructionStrategy;

        ReconstructionStrategyDefinition(ReconstructionStrategy instance) {
            this.reconstructionStrategy = instance;
        }

        public ReconstructionStrategy getReconstructionStrategy() {
            return reconstructionStrategy;
        }
    }

    static class BestByAverage implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease).reversed());
            return contactStructuralInformation.subList(0, numberOfContacts)
                    .stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "highest";
        }
    }

    static class IgnoreBest implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                                 List<ContactStructuralInformation> contactStructuralInformation,
                                                                 int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease).reversed());
            List<ContactStructuralInformation> bestContacts = contactStructuralInformation.subList(0, numberOfContacts);
            Collections.shuffle(contactStructuralInformation);
            return contactStructuralInformation.stream()
                    .filter(contact -> !bestContacts.contains(contact))
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "ignore highest";
        }
    }

    static class Random implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                                 List<ContactStructuralInformation> contactStructuralInformation,
                                                                 int numberOfContacts) {
            Collections.shuffle(contactStructuralInformation);
            return contactStructuralInformation.subList(0, numberOfContacts)
                    .stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "random";
        }
    }

    static class WorstByAverage implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain, List<ContactStructuralInformation> contactStructuralInformation, int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease));
            return contactStructuralInformation.subList(0, numberOfContacts)
                    .stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "lowest";
        }
    }

    /**
     * Select best contacts and introduce 50% non-native contacts additionally.
     */
    static class NonNative implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                                 List<ContactStructuralInformation> contactStructuralInformation,
                                                                 int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease).reversed());
            List<ContactStructuralInformation> information = contactStructuralInformation.subList(0, numberOfContacts);
            List<Pair<Integer, Integer>> contacts = information.stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());

            List<AminoAcid> aminoAcids = chain.getAminoAcids();
            List<Pair<Integer, Integer>> nonNativeContacts = new ArrayList<>();
            // additional non-native contacts
            int numberOfAdditionalContacts = (int) (0.5 * numberOfContacts);
            for (int i = 0; i < numberOfAdditionalContacts; i++) {
                Pair<Integer, Integer> nonNativeContact = null;
                while (nonNativeContact == null) {
                    Collections.shuffle(aminoAcids);
                    AminoAcid aminoAcid1 = aminoAcids.get(0);
                    AminoAcid aminoAcid2 = aminoAcids.get(1);
                    Pair<Integer, Integer> potentialContact = new Pair<>(aminoAcid1.getResidueIdentifier().getResidueNumber(),
                            aminoAcid2.getResidueIdentifier().getResidueNumber());
                    // real contact with enough sequence separation
                    if ((Math.abs(potentialContact.getLeft() - potentialContact.getRight()) > 5) &&
                            // not in native contacts present
                            (!contacts.contains(potentialContact) && !contacts.contains(potentialContact.flip())) &&
                            // not in non-native contacts present
                            (!nonNativeContacts.contains(potentialContact) && !nonNativeContacts.contains(potentialContact.flip()))) {
                        nonNativeContact = potentialContact;
                    }
                }
                nonNativeContacts.add(nonNativeContact);
            }

            return Stream.concat(contacts.stream(), nonNativeContacts.stream())
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "non-native";
        }
    }
}

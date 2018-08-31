package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.align.impl.TMAlignService;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.efr.model.ContactDistanceBin;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.feature.interaction.PLIPInteraction;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinition;
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
    private static final double DEFAULT_COVERAGE = 0.3;
    private static final int REDUNDANCY = 10;
    private static final Path OUTPUT_PATH = Paths.get("/home/sb/strategy-detailed.csv");
    private static final TMAlignService TM_ALIGN_SERVICE = TMAlignService.getInstance();
    private static final ContactDefinition contactDefinition = ContactDefinitionFactory.createAlphaCarbonContactDefinition(8.0);

    private static FileWriter fileWriter;
    private static ExecutorService executorService;

    public static void main(String[] args) throws IOException {
        fileWriter = new FileWriter(OUTPUT_PATH.toFile());
        fileWriter.write("id,strategy,rmsd" + System.lineSeparator());
        executorService = Executors.newFixedThreadPool(16);

        TestUtils.getResourceAsStream("data/efr.list")
                .map(line -> line.split(";"))
                .map(ExplorerChain::new)
                .forEach(A03_ReconstructByVariousStrategy::handleChain);
    }

    private static void handleChain(ExplorerChain explorerChain) {
        logger.info("[{}] starting job",
                explorerChain.getStfId());

        try {
            Chain nativeChain = explorerChain.getChain();
            Path nativeChainPath = Files.createTempFile("nativechain-", ".pdb");
            Files.write(nativeChainPath, nativeChain.getPdbRepresentation().getBytes());

            List<ContactStructuralInformation> contactStructuralInformation = explorerChain.getContacts();

            // annotate with PLIP data
            PLIPInteractionContainer plipInteractionContainer = nativeChain.getFeature(PLIPInteractionContainer.class);
            for(ContactStructuralInformation csi : contactStructuralInformation) {
                AminoAcid aminoAcid1 = nativeChain.select()
                        .residueNumber(csi.getResidueIdentifier1())
                        .asAminoAcid();
                AminoAcid aminoAcid2 = nativeChain.select()
                        .residueNumber(csi.getResidueIdentifier2())
                        .asAminoAcid();
                if(plipInteractionContainer.getHydrogenBonds()
                        .stream()
                        .anyMatch(hydrogenBond -> isContact(hydrogenBond, aminoAcid1, aminoAcid2))) {
                    csi.markAsHydrogenBond();
                }
                if(plipInteractionContainer.getHydrophobicInteractions()
                                .stream()
                                .anyMatch(hydrophobicInteraction -> isContact(hydrophobicInteraction, aminoAcid1, aminoAcid2))) {
                    csi.markAsHydrophobicInteraction();
                }
            }

            int numberOfNativeContacts = contactStructuralInformation.size();
            int numberOfContactsToSelect = (int) (numberOfNativeContacts * DEFAULT_COVERAGE);

            List<ReconstructionContactMap> contactMaps = Stream.of(ReconstructionStrategyDefinition.values())
                    .map(ReconstructionStrategyDefinition::getReconstructionStrategy)
                    .flatMap(reconstructionStrategy -> IntStream.range(0, REDUNDANCY)
                            .mapToObj(i -> {
                                ReconstructionContactMap contactMap = reconstructionStrategy.composeReconstructionContactMap(nativeChain,
                                        contactStructuralInformation,
                                        numberOfContactsToSelect);

                                contactMap.setName(reconstructionStrategy.getName() + "-" + (i + 1));

                                return contactMap;
                            }))
                    .collect(Collectors.toList());

            Map<String, List<Future<ReconstructionResult>>> reconstructionFutures = new HashMap<>();
            for (ReconstructionContactMap contactMap : contactMaps) {
                String name = contactMap.getName().split("-")[0];
                logger.info("[{}] handling contact map definition {}",
                        explorerChain.getStfId(),
                        name);

                if(!reconstructionFutures.containsKey(name)) {
                    reconstructionFutures.put(name, new ArrayList<>());
                }

                List<Future<ReconstructionResult>> bin = reconstructionFutures.get(name);

                bin.add(executorService.submit(new ConfoldServiceWorker("/home/sb/programs/confold_v1.0/confold.pl",
                        contactMap.getSequence(),
                        contactMap.getSecondaryStructureElements(),
                        contactMap.getCaspRRRepresentation(),
                        contactMap.getConfoldRRType())));
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
                    logger.info("[{}][{}] {} reconstructs in bin",
                            explorerChain.getStfId(),
                            name,
                            reconstructions.size());
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

                    logger.info("[{}][{}] {} alignments in bin",
                            explorerChain.getStfId(),
                            name,
                            alignmentResults.size());

                    if (alignmentResults.isEmpty()) {
                        throw new ComputationException("tmalign did not yield any alignments");
                    }

                    for(TMAlignAlignmentResult alignmentResult : alignmentResults) {
                        double rmsd = alignmentResult.getRootMeanSquareDeviation().getScore();
                        String line = explorerChain.getStfId() + "," + name + "," + rmsd;
                        logger.info("[{}][{}] {}",
                                explorerChain.getStfId(),
                                name,
                                line);
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

    private static boolean isContact(PLIPInteraction interaction, AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return (interaction.getPartner1().equals(aminoAcid1) && interaction.getPartner2().equals(aminoAcid2)) ||
                (interaction.getPartner2().equals(aminoAcid1) && interaction.getPartner1().equals(aminoAcid2));
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
                            .collect(Collectors.toList()),
                    contactDefinition);
        }

        String getName();
    }

    enum ReconstructionStrategyDefinition {
        // initial strategies
//        RANDOM(new Random()),
//        BEST_BY_AVERAGE(new BestByAverage()),
//        WORST_BY_AVERAGE(new WorstByAverage()),
//        NON_NATIVE(new BestNativeNonNativeSplit(0, 100)),
//        BEST25_NON_NATIVE75(new BestNativeNonNativeSplit(25, 75)),
//        BEST50_NON_NATIVE50(new BestNativeNonNativeSplit(50, 50)),
//        BEST75_NON_NATIVE25(new BestNativeNonNativeSplit(75, 25)),
//        WORST25_NON_NATIVE75(new WorstNativeNonNativeSplit(25, 75)),
//        WORST50_NON_NATIVE50(new WorstNativeNonNativeSplit(50, 50)),
//        WORST75_NON_NATIVE25(new WorstNativeNonNativeSplit(75, 25)),
//        SHORT(new Short()),
//        LONG(new Long()),
//        HYDROGEN(new HydrogenBond()),
//        HYDROPHOBIC(new HydrophobicInteraction()),
//        // more fine-grained assessment of FP
//        BEST55_NON_NATIVE45(new BestNativeNonNativeSplit(55, 45)),
//        BEST60_NON_NATIVE40(new BestNativeNonNativeSplit(60, 40)),
//        BEST65_NON_NATIVE35(new BestNativeNonNativeSplit(65, 35)),
//        BEST70_NON_NATIVE30(new BestNativeNonNativeSplit(70, 30)),
//        BEST80_NON_NATIVE20(new BestNativeNonNativeSplit(80, 20)),
//        BEST85_NON_NATIVE15(new BestNativeNonNativeSplit(85, 15)),
//        BEST90_NON_NATIVE10(new BestNativeNonNativeSplit(90, 10)),
//        BEST95_NON_NATIVE5(new BestNativeNonNativeSplit(95, 5)),
//        WORST55_NON_NATIVE45(new WorstNativeNonNativeSplit(55, 45)),
//        WORST60_NON_NATIVE40(new WorstNativeNonNativeSplit(60, 40)),
//        WORST65_NON_NATIVE35(new WorstNativeNonNativeSplit(65, 35)),
//        WORST70_NON_NATIVE30(new WorstNativeNonNativeSplit(70, 30)),
//        WORST80_NON_NATIVE20(new WorstNativeNonNativeSplit(80, 20)),
//        WORST85_NON_NATIVE15(new WorstNativeNonNativeSplit(85, 15)),
//        WORST90_NON_NATIVE10(new WorstNativeNonNativeSplit(90, 10)),
//        WORST95_NON_NATIVE5(new WorstNativeNonNativeSplit(95, 5)),
//        RANDOM25_NON_NATIVE75(new RandomNativeNonNativeSplit(25, 75)),
//        RANDOM50_NON_NATIVE50(new RandomNativeNonNativeSplit(50, 50)),
//        RANDOM55_NON_NATIVE45(new RandomNativeNonNativeSplit(55, 45)),
//        RANDOM60_NON_NATIVE40(new RandomNativeNonNativeSplit(60, 40)),
//        RANDOM65_NON_NATIVE35(new RandomNativeNonNativeSplit(65, 35)),
//        RANDOM70_NON_NATIVE30(new RandomNativeNonNativeSplit(70, 30)),
//        RANDOM75_NON_NATIVE30(new RandomNativeNonNativeSplit(75, 25)),
//        RANDOM80_NON_NATIVE20(new RandomNativeNonNativeSplit(80, 20)),
//        RANDOM85_NON_NATIVE15(new RandomNativeNonNativeSplit(85, 15)),
//        RANDOM90_NON_NATIVE10(new RandomNativeNonNativeSplit(90, 10)),
//        RANDOM95_NON_NATIVE5(new RandomNativeNonNativeSplit(95, 5)),
//        // last bin which is quite uninformative
//        BEST45_NON_NATIVE55(new BestNativeNonNativeSplit(45, 55)),
//        BEST40_NON_NATIVE60(new BestNativeNonNativeSplit(40, 60)),
//        BEST35_NON_NATIVE65(new BestNativeNonNativeSplit(35, 65)),
//        BEST30_NON_NATIVE70(new BestNativeNonNativeSplit(30, 70)),
//        BEST20_NON_NATIVE80(new BestNativeNonNativeSplit(20, 80)),
//        BEST15_NON_NATIVE85(new BestNativeNonNativeSplit(15, 85)),
//        BEST10_NON_NATIVE90(new BestNativeNonNativeSplit(10, 90)),
//        BEST5_NON_NATIVE95(new BestNativeNonNativeSplit(5, 95)),
//        WORST45_NON_NATIVE55(new WorstNativeNonNativeSplit(45, 55)),
//        WORST40_NON_NATIVE60(new WorstNativeNonNativeSplit(40, 60)),
//        WORST35_NON_NATIVE65(new WorstNativeNonNativeSplit(35, 65)),
//        WORST30_NON_NATIVE70(new WorstNativeNonNativeSplit(30, 70)),
//        WORST20_NON_NATIVE80(new WorstNativeNonNativeSplit(20, 80)),
//        WORST15_NON_NATIVE85(new WorstNativeNonNativeSplit(15, 85)),
//        WORST10_NON_NATIVE90(new WorstNativeNonNativeSplit(10, 90)),
//        WORST5_NON_NATIVE95(new WorstNativeNonNativeSplit(5, 95)),
//        RANDOM45_NON_NATIVE55(new RandomNativeNonNativeSplit(45, 55)),
//        RANDOM40_NON_NATIVE60(new RandomNativeNonNativeSplit(40, 60)),
//        RANDOM35_NON_NATIVE65(new RandomNativeNonNativeSplit(35, 65)),
//        RANDOM30_NON_NATIVE70(new RandomNativeNonNativeSplit(30, 70)),
//        RANDOM20_NON_NATIVE80(new RandomNativeNonNativeSplit(20, 80)),
//        RANDOM15_NON_NATIVE85(new RandomNativeNonNativeSplit(15, 85)),
//        RANDOM10_NON_NATIVE90(new RandomNativeNonNativeSplit(10, 90)),
//        RANDOM5_NON_NATIVE95(new RandomNativeNonNativeSplit(5, 95)),
        // fine-grained bins up to 10%
        BEST99_NON_NATIVE1(new BestNativeNonNativeSplit(99, 1)),
        BEST98_NON_NATIVE2(new BestNativeNonNativeSplit(98, 2)),
        BEST97_NON_NATIVE3(new BestNativeNonNativeSplit(97, 3)),
        BEST96_NON_NATIVE4(new BestNativeNonNativeSplit(96, 4)),
        BEST94_NON_NATIVE6(new BestNativeNonNativeSplit(94, 6)),
        BEST93_NON_NATIVE7(new BestNativeNonNativeSplit(93, 7)),
        BEST92_NON_NATIVE8(new BestNativeNonNativeSplit(92, 8)),
        BEST91_NON_NATIVE9(new BestNativeNonNativeSplit(91, 9)),
        RANDOM99_NON_NATIVE1(new RandomNativeNonNativeSplit(99, 1)),
        RANDOM98_NON_NATIVE2(new RandomNativeNonNativeSplit(98, 2)),
        RANDOM97_NON_NATIVE3(new RandomNativeNonNativeSplit(97, 3)),
        RANDOM96_NON_NATIVE4(new RandomNativeNonNativeSplit(96, 4)),
        RANDOM94_NON_NATIVE6(new RandomNativeNonNativeSplit(94, 6)),
        RANDOM93_NON_NATIVE7(new RandomNativeNonNativeSplit(93, 7)),
        RANDOM92_NON_NATIVE8(new RandomNativeNonNativeSplit(92, 8)),
        RANDOM91_NON_NATIVE9(new RandomNativeNonNativeSplit(91, 9)),
        ;

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
            return contactStructuralInformation.stream()
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "best";
        }
    }

    static class Random implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            Collections.shuffle(contactStructuralInformation);
            List<Pair<Integer, Integer>> pairs = contactStructuralInformation.stream()
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
            logger.info("[{}] selected {} random contacts, aim: {}",
                    chain.getChainIdentifier().getProteinIdentifier().getPdbId(),
                    pairs.size(),
                    numberOfContacts);
            return pairs;
        }

        @Override
        public String getName() {
            return "random";
        }
    }

    static class WorstByAverage implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease));
            return contactStructuralInformation.subList(0, numberOfContacts)
                    .stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "worst";
        }
    }

    /**
     * Select some percentage best contacts and introduce some percentage non-native contacts additionally.
     */
    static class BestNativeNonNativeSplit implements ReconstructionStrategy {
        private final int nativePercentage;
        private final int nonNativePercentage;

        BestNativeNonNativeSplit(int nativePercentage, int nonNativePercentage) {
            this.nativePercentage = nativePercentage;
            this.nonNativePercentage = nonNativePercentage;
        }

        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease).reversed());
            List<Pair<Integer, Integer>> nativeContacts = contactStructuralInformation.stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());

            List<Pair<Integer, Integer>> selectedNativeContacts = nativeContacts.stream()
                    .limit((long) (nativePercentage * 0.01 * numberOfContacts))
                    .collect(Collectors.toList());

            List<AminoAcid> aminoAcids = chain.getAminoAcids();
            List<Pair<Integer, Integer>> nonNativeContacts = new ArrayList<>();
            // additional non-native contacts
            int numberOfAdditionalContacts = (int) (nonNativePercentage * 0.01 * numberOfContacts);
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
                            (!nativeContacts.contains(potentialContact) && !nativeContacts.contains(potentialContact.flip())) &&
                            // not in non-native contacts present
                            (!nonNativeContacts.contains(potentialContact) && !nonNativeContacts.contains(potentialContact.flip()))) {
                        nonNativeContact = potentialContact;
                    }
                }
                nonNativeContacts.add(nonNativeContact);
            }

            return Stream.concat(selectedNativeContacts.stream(), nonNativeContacts.stream())
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "best" + nativePercentage + "_nonnative" + nonNativePercentage;
        }
    }

    static class WorstNativeNonNativeSplit implements ReconstructionStrategy {
        private final int nativePercentage;
        private final int nonNativePercentage;

        WorstNativeNonNativeSplit(int nativePercentage, int nonNativePercentage) {
            this.nativePercentage = nativePercentage;
            this.nonNativePercentage = nonNativePercentage;
        }

        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            contactStructuralInformation.sort(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease));
            List<Pair<Integer, Integer>> nativeContacts = contactStructuralInformation.stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());

            List<Pair<Integer, Integer>> selectedNativeContacts = nativeContacts.stream()
                    .limit((long) (nativePercentage * 0.01 * numberOfContacts))
                    .collect(Collectors.toList());

            List<AminoAcid> aminoAcids = chain.getAminoAcids();
            List<Pair<Integer, Integer>> nonNativeContacts = new ArrayList<>();
            // additional non-native contacts
            int numberOfAdditionalContacts = (int) (nonNativePercentage * 0.01 * numberOfContacts);
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
                            (!nativeContacts.contains(potentialContact) && !nativeContacts.contains(potentialContact.flip())) &&
                            // not in non-native contacts present
                            (!nonNativeContacts.contains(potentialContact) && !nonNativeContacts.contains(potentialContact.flip()))) {
                        nonNativeContact = potentialContact;
                    }
                }
                nonNativeContacts.add(nonNativeContact);
            }

            return Stream.concat(selectedNativeContacts.stream(), nonNativeContacts.stream())
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "worst" + nativePercentage + "_nonnative" + nonNativePercentage;
        }
    }

    static class RandomNativeNonNativeSplit implements ReconstructionStrategy {
        private final int nativePercentage;
        private final int nonNativePercentage;

        RandomNativeNonNativeSplit(int nativePercentage, int nonNativePercentage) {
            this.nativePercentage = nativePercentage;
            this.nonNativePercentage = nonNativePercentage;
        }

        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            Collections.shuffle(contactStructuralInformation);
            List<Pair<Integer, Integer>> nativeContacts = contactStructuralInformation.stream()
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());

            List<Pair<Integer, Integer>> selectedNativeContacts = nativeContacts.stream()
                    .limit((long) (nativePercentage * 0.01 * numberOfContacts))
                    .collect(Collectors.toList());

            List<AminoAcid> aminoAcids = chain.getAminoAcids();
            List<Pair<Integer, Integer>> nonNativeContacts = new ArrayList<>();
            // additional non-native contacts
            int numberOfAdditionalContacts = (int) (nonNativePercentage * 0.01 * numberOfContacts);
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
                            (!nativeContacts.contains(potentialContact) && !nativeContacts.contains(potentialContact.flip())) &&
                            // not in non-native contacts present
                            (!nonNativeContacts.contains(potentialContact) && !nonNativeContacts.contains(potentialContact.flip()))) {
                        nonNativeContact = potentialContact;
                    }
                }
                nonNativeContacts.add(nonNativeContact);
            }

            return Stream.concat(selectedNativeContacts.stream(), nonNativeContacts.stream())
                    .collect(Collectors.toList());
        }

        @Override
        public String getName() {
            return "random" + nativePercentage + "_nonnative" + nonNativePercentage;
        }
    }

    static class Short implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            List<Pair<Integer, Integer>> pairs = contactStructuralInformation.stream()
                    .filter(csi -> csi.getContactDistanceBin() == ContactDistanceBin.SHORT)
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
            logger.info("[{}] selected {} short contacts, target: {}",
                    chain.getChainIdentifier().getProteinIdentifier().getPdbId(),
                    pairs.size(),
                    numberOfContacts);
            return pairs;
        }

        @Override
        public String getName() {
            return "short";
        }
    }

    static class Long implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            List<Pair<Integer, Integer>> pairs = contactStructuralInformation.stream()
                    .filter(csi -> csi.getContactDistanceBin() == ContactDistanceBin.LONG)
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
            logger.info("[{}] selected {} long contacts, target: {}",
                    chain.getChainIdentifier().getProteinIdentifier().getPdbId(),
                    pairs.size(),
                    numberOfContacts);
            return pairs;
        }

        @Override
        public String getName() {
            return "long";
        }
    }

    static class HydrogenBond implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            List<Pair<Integer, Integer>> pairs = contactStructuralInformation.stream()
                    .filter(ContactStructuralInformation::isHydrogenBond)
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
            logger.info("[{}] selected {} hydrogen bonds, target: {}",
                    chain.getChainIdentifier().getProteinIdentifier().getPdbId(),
                    pairs.size(),
                    numberOfContacts);
            return pairs;
        }

        @Override
        public String getName() {
            return "hydrogen";
        }
    }

    static class HydrophobicInteraction implements ReconstructionStrategy {
        @Override
        public List<Pair<Integer, Integer>> selectContacts(Chain chain,
                                                           List<ContactStructuralInformation> contactStructuralInformation,
                                                           int numberOfContacts) {
            List<Pair<Integer, Integer>> pairs = contactStructuralInformation.stream()
                    .filter(ContactStructuralInformation::isHydrophobicInteraction)
                    .limit(numberOfContacts)
                    .map(contact -> new Pair<>(contact.getResidueIdentifier1(), contact.getResidueIdentifier2()))
                    .collect(Collectors.toList());
            logger.info("[{}] selected {} hydrophobic interactions, target: {}",
                    chain.getChainIdentifier().getProteinIdentifier().getPdbId(),
                    pairs.size(),
                    numberOfContacts);
            return pairs;
        }

        @Override
        public String getName() {
            return "hydrophobic";
        }
    }
}

package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.model.BaselineReconstruction;
import de.bioforscher.jstructure.si.model.ContactTogglingReconstruction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class StructuralInformationService {
    private static final Logger logger = LoggerFactory.getLogger(StructuralInformationService.class);
    private static final StructuralInformationService INSTANCE = new StructuralInformationService();
    /**
     * How many contacts to select for each iteration.
     */
    //TODO could sample some frequencies to normalize structural change expected from one contact
    private static final double BASELINE_FREQUENCY = 0.3;
    /**
     * How many iterations to perform.
     */
    private static final int COVERAGE = 10;

    private StructuralInformationService() {
        logger.info("{} setup up to sample {}% of contacts with coverage {}",
                getClass().getSimpleName(),
                StandardFormat.formatToInteger(BASELINE_FREQUENCY * 100),
                COVERAGE);
    }

    public static StructuralInformationService getInstance() {
        return INSTANCE;
    }

    public void process(Chain referenceChain,
                        String jobName,
                        Path confoldPath,
                        Path tmalignPath,
                        Path outputDirectory,
                        int numberOfThreads) {
        try {
            ExecutorService executorService = Executors.newFixedThreadPool(numberOfThreads);
            List<AminoAcid> aminoAcids = referenceChain.aminoAcids().collect(Collectors.toList());
            int numberOfAminoAcids = aminoAcids.size();
            ReconstructionContactMap fullMap = ReconstructionContactMap.createReconstructionContactMap(referenceChain);
            List<Pair<AminoAcid, AminoAcid>> contacts = fullMap.getLongRangeContacts();
            int numberOfLongRangeContacts = contacts.size();

            String sequence = fullMap.getSequence();
            String secondaryStructure = fullMap.getSecondaryStructureElements();

            // report general statistics
            logger.info("evaluating structural information of chain with {} residues and {} long-range contacts",
                    numberOfAminoAcids,
                    numberOfLongRangeContacts);
            logger.info("present long-range contacts:{}{}",
                    System.lineSeparator(),
                    composeConsoleOutput(fullMap.getLongRangeContacts()));
            logger.info("coverage per position:{}{}",
                    System.lineSeparator(),
                    aminoAcids.stream()
                            .map(fullMap::getNonLocalContactsOf)
                            .mapToInt(List::size)
                            .mapToObj(size -> {
                                if (size < 10) {
                                    return String.valueOf(size);
                                } else {
                                    return "*";
                                }
                            })
                            .collect(Collectors.joining()));


            // write reference structure
            Path referenceChainStructurePath = Files.createTempFile("confoldservice-ref", ".pdb");
            Files.write(referenceChainStructurePath, referenceChain.getPdbRepresentation().getBytes());

            // create sampling containers
            List<Future<BaselineReconstruction>> baselineFutures = new ArrayList<>();
            for (int iteration = 0; iteration < COVERAGE; iteration++) {
                baselineFutures.add(executorService.submit(new BaselineReconstruction(iteration,
                        referenceChainStructurePath,
                        referenceChain,
                        fullMap,
                        sequence,
                        secondaryStructure,
                        BASELINE_FREQUENCY,
                        confoldPath.toFile().getAbsolutePath(),
                        tmalignPath.toFile().getAbsolutePath(),
                        outputDirectory,
                        jobName)));
            }

            List<BaselineReconstruction> baselineReconstructions = baselineFutures.stream()
                    .map(this::getBaselineFuture)
                    .collect(Collectors.toList());

            // create toggling reconstructions
            int numberOfCombinations = baselineReconstructions.size() * contacts.size();
            logger.info("reconstructing individual maps by contact toggling - evaluating {} combinations",
                    numberOfCombinations);
            int counter = 1;
            List<Future<ContactTogglingReconstruction>> contactToggleFutures = new ArrayList<>();
            for (Pair<AminoAcid, AminoAcid> contact : contacts) {
                for (BaselineReconstruction baselineReconstruction : baselineReconstructions) {
                    contactToggleFutures.add(executorService.submit(baselineReconstruction.createContactTogglingReconstruction(contact,
                            counter,
                            numberOfCombinations)));
                    counter++;
                }
            }

            Path outputPath = outputDirectory.resolve(jobName + ".out");
            FileWriter outputWriter = new FileWriter(outputPath.toFile());
            contactToggleFutures.stream()
                    .map(this::getContactToggleFuture)
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .map(contactTogglingReconstruction -> contactTogglingReconstruction.getContactToToggle() + "\t" +
                            contactTogglingReconstruction.isContactWasRemoved() + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getBaselineReconstruction().getAverageRmsd()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getBaselineReconstruction().getAverageTmScore()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getBaselineReconstruction().getAverageQ()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getAverageRmsd()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getAverageTmScore()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getAverageQ()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getDecreaseRmsd()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getIncreaseTMScore()) + "\t" +
                            StandardFormat.format(contactTogglingReconstruction.getIncreaseQ()) + System.lineSeparator())
                    .forEach(line -> {
                        try {
                            outputWriter.write(line);
                            outputWriter.flush();
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    });

            // clean up
            outputWriter.close();
            Files.delete(referenceChainStructurePath);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private Optional<ContactTogglingReconstruction> getContactToggleFuture(Future<ContactTogglingReconstruction> future) {
        try {
            return Optional.of(future.get());
        } catch (Exception e) {
            logger.warn("computation of toggling results failed",
                    e);
            return Optional.empty();
        }
    }

    private BaselineReconstruction getBaselineFuture(Future<BaselineReconstruction> future) {
        try {
            return future.get();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private String composeConsoleOutput(List<Pair<AminoAcid, AminoAcid>> contacts) {
        return contacts.stream()
                .map(contact -> contact.getLeft() + " <-> " + contact.getRight())
                .collect(Collectors.joining(System.lineSeparator()));
    }
}

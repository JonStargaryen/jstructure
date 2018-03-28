package de.bioforscher.jstructure.si.model;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.align.impl.TMAlignService;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.align.result.score.RootMeanSquareDeviation;
import de.bioforscher.jstructure.align.result.score.TemplateModelingScore;
import de.bioforscher.jstructure.graph.ReconstructionContactMap;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.ConfoldServiceWorker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

public class BaselineReconstruction implements Callable<BaselineReconstruction> {
    private static final Logger logger = LoggerFactory.getLogger(BaselineReconstruction.class);
    private static final TMAlignService TM_ALIGN_SERVICE = TMAlignService.getInstance();

    private final int iteration;
    private final Path referenceChainPath;
    private final Chain referenceChain;
    private final ReconstructionContactMap fullMap;
    private final String sequence;
    private final String secondaryStructure;
    private final double baselineFrequency;
    private final String serviceLocation;

    private ReconstructionContactMap sampledMap;
    private List<Chain> sampledReconstructions;
    private double averageRmsd;
    private double averageTmScore;
    private double averageQ;

    public BaselineReconstruction(int iteration,
                                  Path referenceChainPath,
                                  Chain referenceChain,
                                  ReconstructionContactMap fullMap,
                                  String sequence,
                                  String secondaryStructure,
                                  double baselineFrequency,
                                  String serviceLocation) {
        this.iteration = iteration;
        this.referenceChainPath = referenceChainPath;
        this.referenceChain = referenceChain;
        this.fullMap = fullMap;
        this.sequence = sequence;
        this.secondaryStructure = secondaryStructure;
        this.baselineFrequency = baselineFrequency;
        this.serviceLocation = serviceLocation;
    }

    @Override
    public BaselineReconstruction call() throws Exception {
        logger.info("submitting baseline reconstruction job with id {}",
                iteration);

        // create sampling of full map
        List<Pair<AminoAcid, AminoAcid>> fullContacts = fullMap.getLongRangeContacts();
        Collections.shuffle(fullContacts);
        int numberOfContactsToSelect = (int) (fullContacts.size() * baselineFrequency);
        List<Pair<AminoAcid, AminoAcid>> sampledContacts = fullContacts.subList(0, numberOfContactsToSelect);
        this.sampledMap = new ReconstructionContactMap(referenceChain.aminoAcids().collect(Collectors.toList()), sampledContacts);

        // reconstruct sampled baseline map
        sampledReconstructions = new ConfoldServiceWorker(serviceLocation,
                sequence,
                secondaryStructure,
                sampledMap.getCaspRRRepresentation()).call();

        // score baseline models
        computeBaselinePerformance(sampledReconstructions);

        return this;
    }

    private void computeBaselinePerformance(List<Chain> reconstructions) throws IOException, InterruptedException {
        List<TMAlignAlignmentResult> alignmentResults = new ArrayList<>();
        List<ReconstructionContactMap> reconstructionContactMaps = new ArrayList<>();
        for(Chain reconstructedChain : reconstructions) {
            Path reconstructPath = Files.createTempFile("confoldservice-recon", ".pdb");
            Files.write(reconstructPath, reconstructedChain.getPdbRepresentation().getBytes());
            alignmentResults.add(TM_ALIGN_SERVICE.process(new String[] {
                    "/home/sb/programs/tmalign",
                    referenceChainPath.toFile().getAbsolutePath(),
                    reconstructPath.toFile().getAbsolutePath()
            }));
            reconstructionContactMaps.add(ReconstructionContactMap.createReconstructionContactMap(reconstructedChain));
        }

        averageRmsd = alignmentResults.stream()
                .map(TMAlignAlignmentResult::getRootMeanSquareDeviation)
                .mapToDouble(RootMeanSquareDeviation::getScore)
                .average()
                .orElseThrow(() -> new ComputationException("could not generate reconstructs"));
        averageTmScore = alignmentResults.stream()
                .map(TMAlignAlignmentResult::getTemplateModelingScore1)
                .mapToDouble(TemplateModelingScore::getScore)
                .average()
                .orElseThrow(() -> new ComputationException("could not generate reconstructs"));
        averageQ = reconstructionContactMaps.stream()
                .mapToDouble(reconstructContactMap -> computeQ(fullMap, reconstructContactMap))
                .average()
                .orElseThrow(() -> new ComputationException("could not generate reconstructs"));

        logger.info("baseline reconstruction {} - average RMSD: {}, average TM-score: {}, average Q: {}",
                iteration,
                StandardFormat.format(averageRmsd),
                StandardFormat.format(averageTmScore),
                StandardFormat.format(averageQ));
    }

    static double computeQ(ReconstructionContactMap referenceMap, ReconstructionContactMap reconstructMap) {
        int numberOfNativeContacts = referenceMap.getNumberOfContacts();
        return referenceMap.getLongRangeContacts()
                .stream()
                .filter(referencePair -> determineIfMapContainsContact(reconstructMap, referencePair))
                .count() / (double) numberOfNativeContacts;
    }

    private static boolean determineIfMapContainsContact(ReconstructionContactMap map, Pair<AminoAcid, AminoAcid> contact) {
        return map.getLongRangeContacts()
                .stream()
                .anyMatch(pair -> pairsMatchByResidueIdentifier(pair, contact));
    }

    private static boolean pairsMatchByResidueIdentifier(Pair<AminoAcid, AminoAcid> pair1, Pair<AminoAcid, AminoAcid> pair2) {
        String left1 = pair1.getLeft().getResidueIdentifier().toString();
        String right1 = pair1.getRight().getResidueIdentifier().toString();
        String left2 = pair2.getLeft().getResidueIdentifier().toString();
        String right2 = pair2.getRight().getResidueIdentifier().toString();
        return (left2.equals(left1) && right2.equals(right1)) ||
                (right2.equals(left2) && left2.equals(right2));
    }

    public ContactTogglingReconstruction createContactTogglingReconstruction(Pair<AminoAcid, AminoAcid> contactToToggle,
                                                                             int counter,
                                                                             int numberOfCombinations) {
        List<Pair<AminoAcid, AminoAcid>> contacts = new ArrayList<>(sampledMap.getLongRangeContacts());
        boolean contactWasRemoved = determineIfMapContainsContact(sampledMap, contactToToggle);
        if(contactWasRemoved) {
            contacts.removeIf(pair -> pairsMatchByResidueIdentifier(pair, contactToToggle));
        } else {
            contacts.add(contactToToggle);
        }

        return new ContactTogglingReconstruction(this,
                counter,
                numberOfCombinations,
                contactToToggle,
                contactWasRemoved,
                new ReconstructionContactMap(sampledMap.getAminoAcids(), contacts));
    }

    public int getIteration() {
        return iteration;
    }

    public Path getReferenceChainPath() {
        return referenceChainPath;
    }

    public Chain getReferenceChain() {
        return referenceChain;
    }

    public ReconstructionContactMap getFullMap() {
        return fullMap;
    }

    public String getSequence() {
        return sequence;
    }

    public String getSecondaryStructure() {
        return secondaryStructure;
    }

    public double getBaselineFrequency() {
        return baselineFrequency;
    }

    public String getServiceLocation() {
        return serviceLocation;
    }

    public ReconstructionContactMap getSampledMap() {
        return sampledMap;
    }

    public List<Chain> getSampledReconstructions() {
        return sampledReconstructions;
    }

    public double getAverageRmsd() {
        return averageRmsd;
    }

    public double getAverageTmScore() {
        return averageTmScore;
    }

    public double getAverageQ() {
        return averageQ;
    }
}

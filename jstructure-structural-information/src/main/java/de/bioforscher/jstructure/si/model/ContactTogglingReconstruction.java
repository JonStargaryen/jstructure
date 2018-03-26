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
import java.util.List;
import java.util.concurrent.Callable;

public class ContactTogglingReconstruction implements Callable<ContactTogglingReconstruction> {
    private static final Logger logger = LoggerFactory.getLogger(ContactTogglingReconstruction.class);
    private static final TMAlignService TM_ALIGN_SERVICE = TMAlignService.getInstance();

    private final BaselineReconstruction baselineReconstruction;
    private final int counter;
    private final int numberOfCombinations;
    private final Pair<AminoAcid, AminoAcid> contactToToggle;
    private final boolean contactWasRemoved;
    private final ReconstructionContactMap alternateMap;
    private List<Chain> toggledReconstructions;
    private double averageRmsd;
    private double averageTmScore;
    private double averageQ;
    private double deltaRmsd;
    private double deltaTMScore;
    private double deltaQ;

    public ContactTogglingReconstruction(BaselineReconstruction baselineReconstruction,
                                         int counter,
                                         int numberOfCombinations,
                                         Pair<AminoAcid, AminoAcid> contactToToggle,
                                         boolean contactWasRemoved,
                                         ReconstructionContactMap alternateMap) {
        this.baselineReconstruction = baselineReconstruction;
        this.counter = counter;
        this.numberOfCombinations = numberOfCombinations;
        this.contactToToggle = contactToToggle;
        this.contactWasRemoved = contactWasRemoved;
        this.alternateMap = alternateMap;
    }

    @Override
    public ContactTogglingReconstruction call() throws Exception {
        toggledReconstructions = new ConfoldServiceWorker(baselineReconstruction.getServiceLocation(),
                baselineReconstruction.getSequence(),
                baselineReconstruction.getSecondaryStructure(),
                alternateMap.getCaspRRRepresentation())
                .call();

        // score baseline models
        computePerformance(toggledReconstructions);

        return this;
    }

    private void computePerformance(List<Chain> reconstructions) throws IOException, InterruptedException {
        List<TMAlignAlignmentResult> alignmentResults = new ArrayList<>();
        List<ReconstructionContactMap> reconstructionContactMaps = new ArrayList<>();
        for(Chain reconstructedChain : reconstructions) {
            Path reconstructPath = Files.createTempFile("confoldservice-recon", ".pdb");
            Files.write(reconstructPath, reconstructedChain.getPdbRepresentation().getBytes());
            alignmentResults.add(TM_ALIGN_SERVICE.process(new String[] {
                    "tmalign",
                    baselineReconstruction.getReferenceChainPath().toFile().getAbsolutePath(),
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
                .mapToDouble(reconstructContactMap -> BaselineReconstruction.computeQ(baselineReconstruction.getFullMap(), reconstructContactMap))
                .average()
                .orElseThrow(() -> new ComputationException("could not generate reconstructs"));

        logger.info("[{} / {}]: {} reconstruction of contact {}",
                counter,
                numberOfCombinations,
                contactWasRemoved ? "removal" : "addition",
                contactToToggle);
        logger.info("[{} / {}]: average RMSD: {}, average TM-score: {}, average Q: {}",
                counter,
                numberOfCombinations,
                StandardFormat.format(averageRmsd),
                StandardFormat.format(averageTmScore),
                StandardFormat.format(averageQ));

        deltaRmsd = baselineReconstruction.getAverageRmsd() - averageRmsd;
        deltaTMScore = baselineReconstruction.getAverageTmScore() - averageTmScore;
        deltaQ = baselineReconstruction.getAverageQ() - averageQ;
        // invert if contact was removed
        if(contactWasRemoved) {
            deltaRmsd = -deltaRmsd;
            deltaTMScore = -deltaTMScore;
            deltaQ = -deltaQ;
        }
        logger.info("[{} / {}]: delta RMSD: {}, delta TM-score: {}, delta Q: {}",
                counter,
                numberOfCombinations,
                StandardFormat.format(deltaRmsd),
                StandardFormat.format(deltaTMScore),
                StandardFormat.format(deltaQ));
    }

    public BaselineReconstruction getBaselineReconstruction() {
        return baselineReconstruction;
    }

    public Pair<AminoAcid, AminoAcid> getContactToToggle() {
        return contactToToggle;
    }

    public boolean isContactWasRemoved() {
        return contactWasRemoved;
    }

    public ReconstructionContactMap getAlternateMap() {
        return alternateMap;
    }

    public List<Chain> getToggledReconstructions() {
        return toggledReconstructions;
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

    public double getDeltaRmsd() {
        return deltaRmsd;
    }

    public double getDeltaTMScore() {
        return deltaTMScore;
    }

    public double getDeltaQ() {
        return deltaQ;
    }
}

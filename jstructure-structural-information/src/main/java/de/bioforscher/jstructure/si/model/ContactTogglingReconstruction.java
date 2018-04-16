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
    private final Pair<String, String> contactToToggle;
    private final boolean contactWasRemoved;
    private final ReconstructionContactMap alternateMap;
    private double averageRmsd;
    private double averageTmScore;
    private double averageQ;
    private double decreaseRmsd;
    private double increaseTMScore;
    private double increaseQ;

    public ContactTogglingReconstruction(BaselineReconstruction baselineReconstruction,
                                         int counter,
                                         int numberOfCombinations,
                                         Pair<AminoAcid, AminoAcid> contactToToggle,
                                         boolean contactWasRemoved,
                                         ReconstructionContactMap alternateMap) {
        this.baselineReconstruction = baselineReconstruction;
        this.counter = counter;
        this.numberOfCombinations = numberOfCombinations;
        this.contactToToggle = new Pair<>(contactToToggle.getLeft().getResidueIdentifier().toString(),
                contactToToggle.getRight().getResidueIdentifier().toString());
        this.contactWasRemoved = contactWasRemoved;
        this.alternateMap = alternateMap;
    }

    @Override
    public ContactTogglingReconstruction call() throws Exception {
        List<Chain> toggledReconstructions = new ConfoldServiceWorker(baselineReconstruction.getConfoldPath(),
                baselineReconstruction.getSequence(),
                baselineReconstruction.getSecondaryStructure(),
                alternateMap.getCaspRRRepresentation())
                .call();

        // score baseline models
        computePerformance(toggledReconstructions);

        return this;
    }

    private void computePerformance(List<Chain> reconstructions) throws IOException {
        List<TMAlignAlignmentResult> alignmentResults = new ArrayList<>();
        List<ReconstructionContactMap> reconstructionContactMaps = new ArrayList<>();
        List<Path> tmpFiles = new ArrayList<>();

        for(Chain reconstructedChain : reconstructions) {
            Path reconstructPath = Files.createTempFile("confoldservice-recon", ".pdb");
            tmpFiles.add(reconstructPath);
            Files.write(reconstructPath, reconstructedChain.getPdbRepresentation().getBytes());
            alignmentResults.add(TM_ALIGN_SERVICE.process(new String[] {
                    baselineReconstruction.getTmalignPath(),
                    baselineReconstruction.getReferenceChainPath().toFile().getAbsolutePath(),
                    reconstructPath.toFile().getAbsolutePath()
            }));
            reconstructionContactMaps.add(ReconstructionContactMap.createReconstructionContactMap(reconstructedChain));
        }

        averageRmsd = alignmentResults.stream()
                .map(TMAlignAlignmentResult::getRootMeanSquareDeviation)
                .mapToDouble(RootMeanSquareDeviation::getScore)
                .average()
                .orElseThrow(() -> new ComputationException("could not generate toggled reconstructs"));
        averageTmScore = alignmentResults.stream()
                .map(TMAlignAlignmentResult::getTemplateModelingScore1)
                .mapToDouble(TemplateModelingScore::getScore)
                .average()
                .orElseThrow(() -> new ComputationException("could not generate toggled reconstructs"));
        averageQ = reconstructionContactMaps.stream()
                .mapToDouble(reconstructContactMap -> BaselineReconstruction.computeQ(baselineReconstruction.getFullMap(), reconstructContactMap))
                .average()
                .orElseThrow(() -> new ComputationException("could not generate toggled reconstructs"));

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

        if(contactWasRemoved) {
            decreaseRmsd = averageRmsd - baselineReconstruction.getAverageRmsd();
            increaseTMScore = baselineReconstruction.getAverageTmScore() - averageTmScore;
            increaseQ = baselineReconstruction.getAverageQ() - averageQ;
        } else {
            decreaseRmsd = baselineReconstruction.getAverageRmsd() - averageRmsd;
            increaseTMScore = averageTmScore - baselineReconstruction.getAverageTmScore();
            increaseQ = averageQ - baselineReconstruction.getAverageQ();
        }

        logger.info("[{} / {}]: decrease RMSD: {}, increase TM-score: {}, increase Q: {}",
                counter,
                numberOfCombinations,
                StandardFormat.format(decreaseRmsd),
                StandardFormat.format(increaseTMScore),
                StandardFormat.format(increaseQ));

        // cleanup
        for(Path tmpFile : tmpFiles) {
            Files.delete(tmpFile);
        }
    }

    public BaselineReconstruction getBaselineReconstruction() {
        return baselineReconstruction;
    }

    public Pair<String, String> getContactToToggle() {
        return contactToToggle;
    }

    public boolean isContactWasRemoved() {
        return contactWasRemoved;
    }

    public ReconstructionContactMap getAlternateMap() {
        return alternateMap;
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

    public double getDecreaseRmsd() {
        return decreaseRmsd;
    }

    public double getIncreaseTMScore() {
        return increaseTMScore;
    }

    public double getIncreaseQ() {
        return increaseQ;
    }
}

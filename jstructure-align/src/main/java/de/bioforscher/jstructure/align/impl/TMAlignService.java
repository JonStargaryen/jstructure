package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.query.StructureAlignmentQuery;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.align.result.score.RootMeanSquareDeviation;
import de.bioforscher.jstructure.align.result.score.TemplateModelingScore;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.service.ExternalLocalService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

public class TMAlignService extends ExternalLocalService {
    private static final Logger logger = LoggerFactory.getLogger(TMAlignService.class);
    private static final TMAlignService INSTANCE = new TMAlignService();
    private static final String DEFAULT_SERVICE_LOCATION = "tmalign";

    private TMAlignService() {
        super(DEFAULT_SERVICE_LOCATION);
    }

    public static TMAlignService getInstance() {
        return INSTANCE;
    }

    public TMAlignAlignmentResult process(String[] arguments) throws IOException, InterruptedException {
        logger.debug("spawning tmalign process with arguments:{}{}",
                System.lineSeparator(),
                arguments);
        ProcessBuilder processBuilder = new ProcessBuilder(arguments);
        Process process = processBuilder.start();

        List<String> outputLines;
        try(BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                outputLines = br.lines()
                    .collect(Collectors.toList());
        }
        List<String> errorLines;
        try(BufferedReader br = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {
            errorLines = br.lines()
                    .collect(Collectors.toList());
        }

        process.waitFor();

        if(outputLines.stream().anyMatch(line -> line.startsWith("Can not open file:"))) {
            throw new ComputationException("error during tmalign execution:" + System.lineSeparator() +
                    outputLines.stream().collect(Collectors.joining(System.lineSeparator())));
        }
        // errors are not reported in error stream
        if(!errorLines.isEmpty()) {
            throw new ComputationException("error during tmalign execution:" + System.lineSeparator() +
                    errorLines.stream().collect(Collectors.joining(System.lineSeparator())));
        }

        int length1 = 0;
        int length2 = 0;
        int alignedLength = 0;
        RootMeanSquareDeviation rootMeanSquareDeviation = null;
        double seqId = 0;
        TemplateModelingScore templateModelingScore1 = null;
        TemplateModelingScore templateModelingScore2 = null;
        for(String outputLine : outputLines) {
            if(outputLine.startsWith("Length of Chain_1")) {
                length1 = Integer.valueOf(outputLine.split(":")[1].trim().split("\\s+")[0]);
            } else if(outputLine.startsWith("Length of Chain_2")) {
                length2 = Integer.valueOf(outputLine.split(":")[1].trim().split("\\s+")[0]);
            } else if(outputLine.startsWith("Aligned length")) {
                String[] split = outputLine.split("=");
                alignedLength = Integer.valueOf(split[1].split(",")[0].trim());
                rootMeanSquareDeviation = new RootMeanSquareDeviation(Double.valueOf(split[2].split(",")[0].trim()));
                seqId = Double.valueOf(split[4].trim());
            } else if(outputLine.startsWith("TM-score")) {
                double tmscore = Double.valueOf(outputLine.split("=")[1].split("\\(")[0].trim());
                TemplateModelingScore templateModelingScore = new TemplateModelingScore(tmscore);
                if(outputLine.contains("Chain_1")) {
                    templateModelingScore1 = templateModelingScore;
                } else {
                    templateModelingScore2 = templateModelingScore;
                }
            }
        }

        return new TMAlignAlignmentResult(length1,
                length2,
                alignedLength,
                rootMeanSquareDeviation,
                seqId,
                templateModelingScore1,
                templateModelingScore2);
    }

    public TMAlignAlignmentResult process(StructureAlignmentQuery structureAlignmentQuery) {
        try {
            Path referencePath = writeStructureToTemporaryFile(structureAlignmentQuery.getReference());
            Path queryPath = writeStructureToTemporaryFile(structureAlignmentQuery.getQuery());

            String[] arguments = new String[] {
                    getServiceLocation(),
                    referencePath.toFile().toString(),
                    queryPath.toFile().toString()
            };

            TMAlignAlignmentResult result = process(arguments);

            Files.delete(referencePath);
            Files.delete(queryPath);

            return result;
        } catch (IOException | InterruptedException e) {
            throw new ComputationException("could not spawn process", e);
        }
    }
}

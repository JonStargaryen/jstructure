package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.query.StructureAlignmentQuery;
import de.bioforscher.jstructure.align.result.TMAlignAlignmentResult;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.service.ExternalLocalService;

import java.io.IOException;
import java.nio.file.Path;

public class TMAlignService extends ExternalLocalService {
    private static final TMAlignService INSTANCE = new TMAlignService();
    private static final String DEFAULT_SERVICE_LOCATION = "tmalign";

    private TMAlignService() {
        super(DEFAULT_SERVICE_LOCATION);
    }

    public static TMAlignService getInstance() {
        return INSTANCE;
    }

    public TMAlignAlignmentResult process(StructureAlignmentQuery structureAlignmentQuery) {
        try {
            Path referencePath = writeStructureToTemporaryFile(structureAlignmentQuery.getReference());
            Path queryPath = writeStructureToTemporaryFile(structureAlignmentQuery.getQuery());
            Path outputPath = createTemporaryOutputFile();

            executeCommandLineCall(getServiceLocation(),
                    referencePath.toFile().toString(),
                    queryPath.toFile().toString(),
                    ">",
                    outputPath.toFile().toString());

            return parseOutputFile(outputPath);
        } catch (IOException | InterruptedException e) {
            throw new ComputationException("could not spawn process", e);
        }
    }

    private TMAlignAlignmentResult parseOutputFile(Path outputPath) {
        return null;
    }
}

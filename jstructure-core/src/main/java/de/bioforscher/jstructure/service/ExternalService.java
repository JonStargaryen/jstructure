package de.bioforscher.jstructure.service;

import de.bioforscher.jstructure.model.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public abstract class ExternalService implements Service {
    private static final Logger logger = LoggerFactory.getLogger(ExternalService.class);
    private final String serviceLocation;
    private final String servicePrefix;

    ExternalService(String serviceLocation) {
        this.serviceLocation = serviceLocation;
        String serviceName = this.getClass().getSimpleName();
        logger.debug("{} is a {} external service with its location setup to '{}'",
                serviceName,
                this instanceof ExternalLocalService ? "local" : "web",
                this.serviceLocation);
        this.servicePrefix = serviceName.toLowerCase();
    }

    @Override
    public String getServiceLocation() {
        return serviceLocation;
    }

    protected Path writeStructureToTemporaryFile(Structure structure) throws IOException {
        Path tmpFile = Files.createTempFile(servicePrefix, ".pdb");
        Files.write(tmpFile, structure.getPdbRepresentation().getBytes());
        return tmpFile;
    }

    protected Path createTemporaryOutputFile() throws IOException {
        return Files.createTempFile(servicePrefix, ".out");
    }

    protected Path createTemporaryOutputDirectory() throws IOException {
        return Files.createTempDirectory(servicePrefix + "-out");
    }
}

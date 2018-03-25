package de.bioforscher.jstructure.service;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public abstract class ExternalLocalService extends ExternalService {
    private static final Logger logger = LoggerFactory.getLogger(ExternalLocalService.class);

    public ExternalLocalService(String serviceLocation) {
        super(validateProposedServiceLocation(serviceLocation));
    }

    private static String validateProposedServiceLocation(String proposedServiceLocation) {
        //TODO check some registry for overwritten service locations
        return proposedServiceLocation;
    }

    protected void executeCommandLineCall(String... arguments) throws IOException, InterruptedException {
        logger.info("spawning process by arguments:{}{}",
                System.lineSeparator(),
                Stream.of(arguments).collect(Collectors.joining(System.lineSeparator())));
        ProcessBuilder processBuilder = new ProcessBuilder(arguments);
        processBuilder.start().waitFor();
    }
}

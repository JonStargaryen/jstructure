package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.service.ExternalLocalService;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;

public class ConfoldServiceWorker extends ExternalLocalService {
    ConfoldServiceWorker(String serviceLocation) {
        super(serviceLocation);
    }

    /**
     * Reconstructs
     * @param sequence this chain's sequence
     * @param secondaryStructure this chain's 3-state SSE annotation (c, H, E) - may be empty or null and will then be ignored
     * @param contacts list of contacts to consider in CASP-RR format
     */
    public Callable<Void> process(String sequence,
                                  String secondaryStructure,
                                  List<String> contacts) {
        try {
            //TODO impl
            executeCommandLineCall(getServiceLocation());
            return null;
        } catch (IOException | InterruptedException e) {
            throw new ComputationException("could not spawn process", e);
        }
    }
}

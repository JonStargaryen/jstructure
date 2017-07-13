package de.bioforscher.jstructure.feature.uniprot;

import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtBlastService;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

/**
 * Bridge to UniProtJAPI.
 * Created by bittrich on 7/13/17.
 */
public class UniProtBridge {
    private static final UniProtBridge INSTANCE = new UniProtBridge();
    private final UniProtService uniProtService;
    private final UniProtBlastService uniProtBlastService;

    private UniProtBridge() {
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        uniProtService = serviceFactoryInstance.getUniProtQueryService();
        uniProtService.start();
        uniProtBlastService = serviceFactoryInstance.getUniProtBlastService();
        uniProtBlastService.start();
    }

    public static UniProtBridge getInstance() {
        return INSTANCE;
    }

    public UniProtBlastService getUniProtBlastService() {
        return uniProtBlastService;
    }

    public UniProtService getUniProtService() {
        return uniProtService;
    }
}

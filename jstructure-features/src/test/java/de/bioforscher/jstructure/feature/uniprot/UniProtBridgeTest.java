package de.bioforscher.jstructure.feature.uniprot;

import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.uniprot.dataservice.client.QueryResult;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;
import uk.ac.ebi.uniprot.dataservice.query.Query;

/**
 * Test integration of the UniProtJAPI.
 * Created by bittrich on 7/13/17.
 */
public class UniProtBridgeTest {
    private UniProtService uniProtService;

    @Before
    public void setup() {
        uniProtService = UniProtBridge.getInstance().getUniProtService();
    }

    @Test
    public void shouldRetrieveSequenceFromUniProtById() throws ServiceException {
        Query query = UniProtQueryBuilder.id("CYC32_DESNO");
        QueryResult<UniProtEntry> result = uniProtService.getEntries(query);
        UniProtEntry entry = result.getFirstResult();
        System.out.println(entry.getSequence());
    }

    @Test
    public void shouldRetrieveSequenceFromUniProtByOldId() throws ServiceException {
        Query query = UniProtQueryBuilder.id("CYC32_DESDN");
        QueryResult<UniProtEntry> result = uniProtService.getEntries(query);
        UniProtEntry entry = result.getFirstResult();
        System.out.println(entry.getSequence());
    }
}
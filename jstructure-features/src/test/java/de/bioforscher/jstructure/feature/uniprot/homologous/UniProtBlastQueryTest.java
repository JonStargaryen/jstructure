package de.bioforscher.jstructure.feature.uniprot.homologous;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtHit;

import java.util.List;

/**
 * Test integration of the UniProtJAPI.
 * Created by bittrich on 7/10/17.
 */
public class UniProtBlastQueryTest {
    @Test
    public void shouldRunBlastQuery() {
        List<UniProtHit> hits = new UniProtHomologyAnnotator().runUniProtBlastService("MES00005665499\n" +
                "MSNHGFAYFFTSYQSLSLDSSSPPPSPHPRAHASSRFPPRARAVASFHTSCKMARTKQTA\n" +
                "RKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRKYQKSTELLIRKLPF\n" +
                "QRLVREIAQDFKTDLRFQSSAVLALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDVQLA\n" +
                "RRIRGERA");
        Assert.assertTrue(hits.size() >= 50);
        hits.stream()
                .map(UniProtHit::getEntry)
                .map(entry -> entry.getPrimaryUniProtAccession() + " " + entry.getFeatures(FeatureType.MUTAGEN).size() +
                        " mutations " + entry.getFeatures(FeatureType.VARIANT).size() + " variants")
                .forEach(System.out::println);
    }
}
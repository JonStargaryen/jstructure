package de.bioforscher.jstructure.feature.interactions;


import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Test;

/**
 * Tests the PLIP parser.
 * Created by bittrich on 2/9/17.
 */
public class PLIPParserTest {
    @Test
    public void shouldAnnotateProtein() {
        String pdbId = "1brr";
        String chainId = "A";
        Chain chain = ProteinParser.source(pdbId).parse()
                .select()
                .chainName(chainId)
                .asChain();
        PLIPParser.parse(chain, PLIPRestServiceQuery.getDocument(pdbId, chainId));
    }
}
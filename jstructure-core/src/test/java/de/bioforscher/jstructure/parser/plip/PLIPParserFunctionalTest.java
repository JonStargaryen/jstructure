package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

/**
 * Tests the PLIP parser.
 * Created by bittrich on 2/9/17.
 */
public class PLIPParserFunctionalTest {
    @Test
    public void shouldAnnotateProtein() {
        String pdbId = "1brr";
        String chainId = "A";
        Chain chain = Selection.on(ProteinParser.source(pdbId).parse())
                .chainName(chainId)
                .asChain();
        PLIPParser.parse(chain, PLIPRestServiceQuery.getPlipResults(pdbId, chainId));
    }
}

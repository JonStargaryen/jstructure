package de.bioforscher.jstructure.feature.interaction;


import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
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
        Chain chain = StructureParser.fromPdbId(pdbId).parse()
                .select()
                .chainName(chainId)
                .asChain();
        PLIPParser.parse(chain, PLIPRestServiceQuery.getIntraMolecularDocument(pdbId, chainId));
    }
}

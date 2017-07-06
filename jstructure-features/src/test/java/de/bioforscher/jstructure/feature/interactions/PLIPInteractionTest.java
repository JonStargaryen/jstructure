package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.testutil.TestUtils;
import org.jsoup.Jsoup;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.stream.Stream;

/**
 * Tests functions of PLIPInteraction instances.
 * Created by bittrich on 6/9/17.
 */
public class PLIPInteractionTest {
    @Test
    public void shouldDistinguishBetweenBackboneAndSideChainInteractions() throws IOException {
        Protein protein = ProteinParser.source("1ar1").parse();
        Chain chain = protein.select()
                .chainName("A")
                .asChain();
        List<PLIPInteraction> interactions = PLIPParser.parse(chain,
                Jsoup.parse(TestUtils.getResourceAsInputStream("plip/1ar1_A.xml"), "UTF-8", ""));
        interactions.forEach(interaction -> {
            System.out.println(interaction);
            // only exactly one can be true
            Boolean a = interaction.isBackboneInteraction();
            Boolean b = interaction.isSideChainInteraction();
            Boolean c = interaction.isMixedInteraction();
            Assert.assertTrue(Stream.of(a, b, c)
                    .filter(Boolean::booleanValue)
                    .count() == 1);
        });
    }
}
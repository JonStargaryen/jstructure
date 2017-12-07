package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Test;

public class EvolutionaryCouplingParserTest {
    @Test
    public void shouldParseHotSpotFile() {
        Structure structure = StructureParser.source("1a64").parse();
        Chain chain = structure.chains().findFirst().get();
        EvolutionaryCouplingParser.parseHotSpotFile(chain, TestUtils.getResourceAsInputStream("couplings/1a64_A_hs.html"));

        chain.aminoAcids()
                .forEach(aminoAcid -> System.out.println(aminoAcid + " > " + aminoAcid.getFeature(EvolutionaryCouplingParser.HotSpotScoring.class).getCumStrength()));
    }
}
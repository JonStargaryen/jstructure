package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;

public class StructuralInformationServiceTest {
    private Chain shortChain;
    private Chain chain1bdd;
    private StructuralInformationService structuralInformationService;

    @Before
    public void setup() {
        shortChain = StructureParser.fromInputStream(TestUtils.getResourceAsInputStream("confold/short.pdb"))
                .parse()
                .getFirstChain();
        chain1bdd = StructureParser.fromPdbId("1bdd").parse().getFirstChain();
        structuralInformationService = StructuralInformationService.getInstance();
    }

    @Test
    @Ignore
    public void shouldCalculateStructuralInformation() throws IOException {
        structuralInformationService.process(shortChain);
    }

    @Test
    @Ignore
    public void shouldCalculateStructuralInformationFor1bdd() throws IOException {
        structuralInformationService.process(chain1bdd);
    }
}
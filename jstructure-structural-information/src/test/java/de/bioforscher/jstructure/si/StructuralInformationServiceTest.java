package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

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
    public void shouldCalculateStructuralInformation() {
        structuralInformationService.process(shortChain, Paths.get("/tmp/short-si.out"));
    }

    @Test
    @Ignore
    public void shouldCalculateStructuralInformationFor1bdd() {
        structuralInformationService.process(chain1bdd, Paths.get("/tmp/1bdd-si.out"));
    }

    @Test
    @Ignore
    public void shouldCalculateStructuralInformationForEarlyFoldingDataset() throws IOException {
        Path outputPath = Paths.get("/home/bittrich/git/phd_sb_repo/data/si/raw/");
        Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/pdb/"))
                .filter(path -> !path.toFile().getName().startsWith("STF0045"))
                .forEach(path -> {
                    Chain chain = StructureParser.fromPath(path)
                            .parse()
                            .getFirstChain();
                    structuralInformationService.process(chain, outputPath.resolve(path.toFile().getName().split("\\.")[0] + ".out"));
                });
    }
}
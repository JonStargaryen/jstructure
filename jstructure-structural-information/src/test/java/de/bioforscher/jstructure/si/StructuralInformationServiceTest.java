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
        structuralInformationService.process(shortChain,
                "shortrun",
                Paths.get("/home/bittrich/programs/confold_v1.0/confold.pl"),
                Paths.get("/home/bittrich/programs/tmalign/tmalign"),
                Paths.get("/tmp/"),
                8,
                0.3);
    }

    @Test
    @Ignore
    public void shouldCalculateStructuralInformationFor1bdd() {
        structuralInformationService.process(chain1bdd,
                "1bddrun",
                Paths.get("/home/bittrich/programs/confold_v1.0/confold.pl"),
                Paths.get("/home/bittrich/programs/tmalign/tmalign"),
                Paths.get("/tmp/"),
                8,
                0.3);
    }

    @Test
    @Ignore
    public void shouldCalculateStructuralInformationForEarlyFoldingDataset() throws IOException {
        Path outputPath = Paths.get("/home/bittrich/si/");
        Files.list(Paths.get("/home/bittrich/si/"))
                .filter(path -> !Files.isDirectory(path))
                .forEach(path -> {
                    Chain chain = StructureParser.fromPath(path)
                            .parse()
                            .getFirstChain();
                    structuralInformationService.process(chain,
                            path.toFile().getName().split("\\.")[0],
                            Paths.get("/home/bittrich/programs/confold_v1.0/confold.pl"),
                            Paths.get("/home/bittrich/programs/tmalign/tmalign"),
                            outputPath,
                            8,
                            0.3);
                });
    }
}
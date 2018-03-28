package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class StructuralInformationRunner {
    private static final StructuralInformationService STRUCTURAL_INFORMATION_SERVICE = StructuralInformationService.getInstance();

    public static void main(String[] args) throws IOException {
        Path workingDirectory = Paths.get(args[0]);
        Files.list(workingDirectory)
                .filter(path -> path.toFile().getName().endsWith(".pdb"))
                .forEach(path -> {
                    Chain chain = StructureParser.fromPath(path)
                            .parse()
                            .getFirstChain();
                    STRUCTURAL_INFORMATION_SERVICE.process(chain, workingDirectory.resolve(path.toFile().getName().split("\\.")[0] + ".out"));
                });
    }
}

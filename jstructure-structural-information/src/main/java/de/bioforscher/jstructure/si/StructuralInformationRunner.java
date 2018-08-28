package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class StructuralInformationRunner {
    private static final Logger logger = LoggerFactory.getLogger(StructuralInformationRunner.class);
    private static final StructuralInformationService STRUCTURAL_INFORMATION_SERVICE = StructuralInformationService.getInstance();

    public static void main(String[] args) throws IOException {
        if(args.length != 4) {
            logger.warn("usage: java -jar distiller.jar [confoldPath] [tmalignPath]  [workingDirectory] [threads]");
            return;
        }

        Path confoldPath = Paths.get(args[0]);
        Path tmalignPath = Paths.get(args[1]);
        Path workingDirectory = Paths.get(args[2]);
        int numberOfThreads = Integer.parseInt(args[3]);
        Files.list(workingDirectory)
                .filter(path -> path.toFile().getName().endsWith(".pdb"))
                .forEach(path -> {
                    logger.info("starting job on {}",
                            path);
                    Chain chain = StructureParser.fromPath(path)
                            .parse()
                            .getFirstChain();
                    STRUCTURAL_INFORMATION_SERVICE.process(chain,
                            path.toFile().getName().split("\\.")[0],
                            confoldPath,
                            tmalignPath,
                            workingDirectory,
                            numberOfThreads);
                });
    }
}

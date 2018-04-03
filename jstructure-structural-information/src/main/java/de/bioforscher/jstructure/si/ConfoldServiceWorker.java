package de.bioforscher.jstructure.si;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

public class ConfoldServiceWorker implements Callable<List<Chain>> {
    private static final Logger logger = LoggerFactory.getLogger(ConfoldServiceWorker.class);
    private final String serviceLocation;
    private final String sequence;
    private final String secondaryStructure;
    private final String contacts;

    public ConfoldServiceWorker(String serviceLocation,
                         String sequence,
                         String secondaryStructure,
                         String contacts) {
        this.serviceLocation = serviceLocation;
        this.sequence = sequence;
        this.secondaryStructure = secondaryStructure;
        this.contacts = contacts;
    }

    @Override
    public List<Chain> call() {
        try {
            Path sequencePath = Files.createTempFile("confoldworker-seq", ".fasta");
            Path mapPath = Files.createTempFile("confoldworker-map", ".rr");
            Path secondaryStructurePath = Files.createTempFile("confoldworker-sse", ".ss");
            Path outputDirectory = Files.createTempDirectory("confoldworker-output");

            Files.write(sequencePath, sequence.getBytes());
            Files.write(mapPath, contacts.getBytes());
            Files.write(secondaryStructurePath, secondaryStructure.getBytes());

            String[] arguments = new String[] {
                    serviceLocation,
                    "-rrtype",
                    "ca",
                    "-seq",
                    sequencePath.toFile().getAbsolutePath(),
                    "-rr",
                    mapPath.toFile().getAbsolutePath(),
                    "-ss",
                    secondaryStructurePath.toFile().getAbsolutePath(),
                    "-o",
                    outputDirectory.toFile().getAbsolutePath()
            };

            logger.debug("spawning confold process with arguments:{}{}",
                    System.lineSeparator(),
                    arguments);
            ProcessBuilder processBuilder = new ProcessBuilder(arguments);
            Process process = processBuilder
//                    .inheritIO()
                    .start();

            process.waitFor();

            List<Chain> chains = Files.list(outputDirectory.resolve("stage1"))
                    .filter(path -> path.toFile().getName().endsWith(".pdb"))
                    .filter(path -> path.toFile().getName().contains("_model"))
                    .map(path -> StructureParser.fromPath(path).parse())
                    .map(Structure::getFirstChain)
                    .collect(Collectors.toList());

            // cleanup
            Files.delete(sequencePath);
            Files.delete(mapPath);
            Files.delete(secondaryStructurePath);
            Files.walk(outputDirectory, FileVisitOption.FOLLOW_LINKS)
                    .sorted(Comparator.reverseOrder())
                    .map(Path::toFile)
                    .forEach(File::delete);

            return chains;
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }
}
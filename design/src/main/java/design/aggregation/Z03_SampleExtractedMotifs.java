package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by S on 29.11.2016.
 */
@Deprecated
public class Z03_SampleExtractedMotifs {
    private static final int limit = 100;

    public static void main(String[] args) throws IOException {
        Stream.of(SequenceMotifDefinition.values())
                .forEach(motif -> {
                    System.out.println(motif);
                    try {
                        Files.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR))
                                .forEach(topologyDir -> {
                                    System.out.println(topologyDir.toFile().getName());
                                    try {
                                        List<Path> motifFiles = Files.list(topologyDir)
                                                .filter(motifFile -> motifFile.toFile().getName().startsWith(motif.name()))
                                                .collect(Collectors.toList());
                                        Collections.shuffle(motifFiles);
                                        motifFiles.stream()
                                                .limit(limit)
                                                .forEach(motifFile -> {
                                                    try {
                                                        Files.copy(motifFile, Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_SAMPLED_DIR + topologyDir.toFile().getName() + "/" + motifFile.toFile().getName()));
                                                    } catch (IOException e) {
                                                        throw new UncheckedIOException(e);
                                                    }
                                                });
                                    } catch (IOException e) {
                                        throw new UncheckedIOException(e);
                                    }
                                });
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}

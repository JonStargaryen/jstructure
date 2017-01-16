package design.visualization;

import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Move all consensus fragments to 1 directory to visualize them by PyMOL.
 * Created by bittrich on 12/15/16.
 */
public class S01_ConsensusFragmentsToPyMOL {
    public static void main(String[] args) throws IOException {
        Map<String, List<Path>> consensusPaths = Files.walk(Paths.get(DesignConstants.FRAGMENT_CLUSTERS_DIR + "trans/"))
                .filter(path -> path.toFile().getName().startsWith(DesignConstants.CLUSTER_CONSENSUS))
                .filter(path -> !path.toFile().getAbsolutePath().contains("/0/"))
                .collect(Collectors.groupingBy(path -> path.toFile().getAbsolutePath().split("/")[6]));

        consensusPaths.entrySet().stream()
                .filter(entry -> entry.getValue().size() > 1)
                .map(Map.Entry::getValue)
                .flatMap(Collection::stream)
                .forEach(path -> {
                    String filepath = path.toFile().getAbsolutePath();
                    System.out.println(filepath);
                    String[] split = filepath.split("/");
                    try {
                        Files.copy(path, Paths.get(DesignConstants.FRAGMENT_CONSENSUS_CLUSTERS + split[5] + "-" + split[6] + "-" + split[7] + DesignConstants.PDB_SUFFIX));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}

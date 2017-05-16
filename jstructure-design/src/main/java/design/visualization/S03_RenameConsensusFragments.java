package design.visualization;

import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Rename to show names conveniently in PyMOL.
 * Created by bittrich on 12/15/16.
 */
public class S03_RenameConsensusFragments {
    public static void main(String[] args) throws IOException {
        Files.list(Paths.get(DesignConstants.FRAGMENT_CONSENSUS_CLUSTERS))
                .filter(path -> path.toFile().getName().startsWith("aligned-"))
                .forEach(path -> {
                    try {
                        Files.copy(path, Paths.get(DesignConstants.FRAGMENT_CONSENSUS_CLUSTERS + "aligned/" +
                                path.toFile().getName().replace("aligned-trans-", "")));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}

package studies.interactions.topology;

import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import studies.StudyConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Are there certain preferences for given interactions to be located in the membrane or outside of it? Or surrounding
 * sequences? Consider if interactions are long-range (i.e. more than 4 positions apart), consider if they are backbone-
 * based or mediated by side-chains.
 *
 * <ul>
 *     <li>load chains</li>
 *     <li>assign interactions</li>
 *     <li>assign topology</li>
 *     <li>assign sequence motifs</li>
 *     <li>write tsv</li>
 * </ul>
 * Created by bittrich on 5/16/17.
 */
public class LinkInteractionsToTopologyAndSequenceMotifs {
    private static final PLIPAnnotator plipAnnotator = new PLIPAnnotator();
    private static final SequenceMotifAnnotator sequenceMotifAnnotator = new SequenceMotifAnnotator();

    public static void main(String[] args) throws IOException {
        Files.list(Paths.get(StudyConstants.GIT + "phd_sb_repo/data/chains/"))
                .forEach(LinkInteractionsToTopologyAndSequenceMotifs::handleFile);
    }

    private static void handleFile(Path path) {
        System.out.println(path);
    }
}

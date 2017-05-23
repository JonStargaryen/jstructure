package studies.interactions.topology;

import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.io.IOException;
import java.nio.file.Files;
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
class LinkInteractionsToTopologyAndSequenceMotifs {
    private static final PLIPAnnotator plipAnnotator = new PLIPAnnotator();
    private static final SequenceMotifAnnotator sequenceMotifAnnotator = new SequenceMotifAnnotator();
    private static final OrientationsOfProteinsInMembranesAnnotator topologyAnnotator = new OrientationsOfProteinsInMembranesAnnotator();

    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(""))
                .limit(5)
                .forEach(LinkInteractionsToTopologyAndSequenceMotifs::handleLine);
    }

    private static void handleLine(String line) {
        System.out.println(line);
        Protein protein = ProteinParser.source(line.split("_")[0]).parse();
        Chain chain = protein.select().chainName(line.split("_")[1]).asChain();

        plipAnnotator.process(protein);
        sequenceMotifAnnotator.process(protein);
        topologyAnnotator.process(protein);
    }
}

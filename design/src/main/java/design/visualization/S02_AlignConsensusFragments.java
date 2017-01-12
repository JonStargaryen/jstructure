package design.visualization;

import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Aligns all consensus fragments relative to one reference structure.
 * Created by bittrich on 12/15/16.
 */
public class S02_AlignConsensusFragments {
    public static void main(String[] args) throws IOException {
        List<Protein> proteins = Files.list(Paths.get(DesignConstants.FRAGMENT_CONSENSUS_CLUSTERS))
                .map(ProteinParser::parsePDBFile)
                .collect(Collectors.toList());

        //TODO manually move first structure
        Protein reference = proteins.remove(0);

        proteins.stream()
                .map(protein -> new SVDSuperimposer().align(reference, protein))
                .forEach(alignmentResult -> {
                    String identifier = alignmentResult.getOriginalQuery().getIdentifier();
                    AtomContainer container = alignmentResult.getOriginalQuery();
                    alignmentResult.transform(container);
                    DesignConstants.write(Paths.get(DesignConstants.FRAGMENT_CONSENSUS_CLUSTERS + "aligned-" + identifier + DesignConstants.PDB_SUFFIX),
                            container.composePDBRecord().getBytes());
                });
    }
}

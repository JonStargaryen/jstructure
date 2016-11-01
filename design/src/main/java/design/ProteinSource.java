package design;

import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.parser.opm.OPMParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Provides access to all proteins of the data set with attached <code>*.opm</code> topology information.
 * Created by S on 31.10.2016.
 */
public class ProteinSource {
    /**
     * The universal function to load protein structures. All proteins of the data set which also provide <code>*.opm</code>
     * files will be returned with topology information attached to them.
     * @return a collection of all
     * @throws IOException when directories or files cannot be found
     */
    public static List<Protein> loadProteins() throws IOException {
        List<Protein> proteins = Files.list(Paths.get(DesignConstants.OPM_RAW_DIR))
                                      .limit(1)
                                      .map(path -> ProteinParser.parsePDBFile(DesignConstants.PDB_DIR +
                                              path.toFile().getName().split("\\.")[0] + DesignConstants.PDB_SUFFIX))
                                      .collect(Collectors.toList());

        proteins.forEach(protein -> {
            System.out.println("parsing information for " + protein.getName());
            OPMParser.parse(protein, Paths.get(DesignConstants.OPM_RAW_DIR + protein.getName() +
                    DesignConstants.OPM_SUFFIX));
            addSecondaryStructureInformation(protein);
            addSequenceMotifInformation(protein);
        });

        return proteins;
    }

    private static void addSequenceMotifInformation(Protein protein) {
        SequenceMotifAnnotator sma = new SequenceMotifAnnotator();
        sma.process(protein);
    }

    private static void addSecondaryStructureInformation(Protein protein) {
        SecondaryStructureAnnotator ssa = new SecondaryStructureAnnotator();
        ssa.process(protein);
    }
}

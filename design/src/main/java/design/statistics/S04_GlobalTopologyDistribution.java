package design.statistics;

import de.bioforscher.jstructure.model.structure.selection.Selection;
import design.parser.opm.OPMParser;
import design.parser.opm.TMHelix;

import java.io.IOException;
import java.util.Map;
import java.util.stream.Collectors;

import static design.DesignConstants.DELIMITER;
import static design.ProteinSource.loadProteins;
import static java.util.Objects.nonNull;

/**
 * Reports the distribution of transmembrane and non-transmembrane regions within the data set.
 * Created by S on 07.12.2016.
 */
@Deprecated
public class S04_GlobalTopologyDistribution {
    public static void main(String[] args) throws IOException {
        // findAny proteins with sequence motif information
        Map<Boolean, Long> count = loadProteins(false, false, true)
                .stream()
                .flatMap(protein -> Selection.on(protein)
                        .aminoAcids()
                        .asFilteredGroups())
                .collect(Collectors.groupingBy(residue -> nonNull(residue.getFeature(TMHelix.class,
                        OPMParser.FeatureNames.TM_HELIX)), Collectors.counting()));

        System.out.println("topology" + DELIMITER + "count");
        System.out.println("transmembrane" + DELIMITER + count.get(true));
        System.out.println("non-transmembrane" + DELIMITER + count.get(false));
    }
}

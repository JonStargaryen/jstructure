package bioforscher.kinks;

import bioforscher.Constants;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Parse KinkFinder results.
 * Created by bittrich on 2/14/17.
 */
public class KinkFinderParser {
    public static final String KINK_FINDER_ANNOTATION = "KINK_FINDER_ANNOTATION";
    private static final String DELIMITER = ",";

    static void parseKinkFinderFile(Protein protein, Path path) {
        List<KinkFinderHelix> helices = Constants.lines(path)
                .filter(line -> !line.startsWith("pdb_code"))
                .map(line -> line.split(DELIMITER))
                .map(KinkFinderHelix::new)
                .collect(Collectors.toList());

        helices.forEach(kinkFinderHelix -> assignHelix(protein, kinkFinderHelix));
        protein.setFeature(KINK_FINDER_ANNOTATION, helices);
    }

    @SuppressWarnings("unchecked")
    private static void assignHelix(Protein protein, KinkFinderHelix kinkFinderHelix) {
        String chainId = kinkFinderHelix.getPdbCode().substring(4, 5);
        IntegerRange range = new IntegerRange(kinkFinderHelix.getHelixStart(), kinkFinderHelix.getHelixEnd());

        Selection.on(protein)
                .chainName(chainId)
                .aminoAcids()
                .residueNumber(range)
                .asFilteredGroups()
                .forEach(group -> {
                    //TODO standardize
                    List<KinkFinderHelix> value = group.getFeature(List.class, KINK_FINDER_ANNOTATION);
                    // entry will be null at first - create list and assign reference
                    if(value == null) {
                        value = new ArrayList<>();
                        group.setFeature(KINK_FINDER_ANNOTATION, value);
                    }
                    value.add(kinkFinderHelix);
                });
    }
}

package de.bioforscher.jstructure.parser.kinkfinder;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Parse KinkFinder results.
 * TODO standardize how parsers behave and how they can be accessed
 * Created by bittrich on 2/14/17.
 */
public class KinkFinderParser {
    public static final String KINK_FINDER_ANNOTATION = "KINK_FINDER_ANNOTATION";
    private static final String DELIMITER = ",";

    public static void parseKinkFinderFile(Protein protein, Path path) throws IOException {
        List<KinkFinderHelix> helices = Files.lines(path)
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
                .forEach(group -> group.addFeatureToList(KINK_FINDER_ANNOTATION, kinkFinderHelix));
    }
}

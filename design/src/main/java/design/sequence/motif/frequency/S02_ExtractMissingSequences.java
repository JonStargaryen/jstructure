package design.sequence.motif.frequency;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

/**
 * Only populated clusters will be extracted by {@link S01_ExtractSequences}. This time, the rarely occurring sequences
 * will be extracted.
 * Created by S on 13.12.2016.
 */
public class S02_ExtractMissingSequences {
    public static void main(String[] args) {
        System.out.println(extractRareSequences("tm", SequenceMotifDefinition.LY6));
    }

    static String extractRareSequences(String topology, SequenceMotifDefinition sequenceMotifDefinition) {
        return DesignConstants.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/"))
                // check that an extracted fragment exists, but the aligned-dir does not contain this entry (as it could
                // not be aligned to any consensus fragment)
                .filter(path -> path.toFile().getName().startsWith(sequenceMotifDefinition.name()))
                .filter(path -> {
                    String fileNameToFind = path.toFile().getName().split("\\.")[0];
                    return DesignConstants.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/"))
                            .map(Path::toFile)
                            .map(File::getName)
                            .noneMatch(fileName -> fileName.startsWith(fileNameToFind));
                })
                .map(path -> ProteinParser.source(path).parse())
                .map(Protein::getAminoAcidSequence)
                .collect(Collectors.joining(System.lineSeparator()));
    }
}

package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import studies.membrane.MembraneConstants;

import java.nio.file.Path;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Aim:
 * Most interactions embedded in sequence motifs are hydrogen bonds. However, are there exceptions?
 *
 * Results:
 * There are 42 water bridges in the data set - seemingly spread evenly among sequence motifs. Nothing particular to
 * gain here.
 *
 * Created by bittrich on 6/22/17.
 */
public class T023_FilterForNonHydrogenBonds {
    public static void main(String[] args) {
        String output = Stream.of(SequenceMotifDefinition.values())
                .flatMap(T023_FilterForNonHydrogenBonds::handleSequenceMotif)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(MembraneConstants.PDBTM_FRAGMENTS_TM_BY_SEQUENCE_MOTIF_PATH.resolve("interactions-non-hydrogen-bond.tsv"),
                output);
    }

    private static Stream<String> handleSequenceMotif(SequenceMotifDefinition sequenceMotifDefinition) {
        Path interactionTsv = MembraneConstants.PDBTM_FRAGMENTS_TM_BY_SEQUENCE_MOTIF_PATH.resolve("interactions-"
                + sequenceMotifDefinition.name() + ".tsv");

        return MembraneConstants.lines(interactionTsv)
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                // filter for non-empty interactions
                .filter(split -> !split[1].equals("-"))
                // filter for non-hydrogen bonds
                .filter(split -> !split[1].equals("Y"))
                .map(split -> sequenceMotifDefinition.name() + "\t" + split[0] + "\t" + split[1] + "\t" + split[2] + "\t" + split[3]);
    }
}

package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.Pair;
import studies.membrane.MembraneConstants;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.Stream;

/**
 * Aim:
 * In the previous step information on interactions within sequence motifs was gathered. Furthermore, for each
 * observation the position of the interaction was stored. E.g. the 1st position interacts with the amino acid at the
 * 4th position via a hydrogen bond.
 * Output is a matrix describing the frequency of interactions for each positional combination. Interactions should
 * always be 'forward' (i.e. O -> N) and the total occurrence count is normalized (i.e. divided by the total number a
 * certain sequence motif occurs in the data set). This still means that the sum of all interactions within a row or
 * column of the matrix will not be 1.0, neither is the upper bound of any entry 1.0 as one sequence motif may embed
 * more than one interaction.
 *
 * Results:
 * Matrices for all sequence motifs. Visualize by R.
 *
 * Created by bittrich on 6/21/17.
 */
public class T022_CondenseFragmentInteractionPositions {
    public static void main(String[] args) {
        Stream.of(SequenceMotifDefinition.values())
                .forEach(T022_CondenseFragmentInteractionPositions::handleSequenceMotif);
    }

    private static void handleSequenceMotif(SequenceMotifDefinition sequenceMotifDefinition) {
        Path interactionTsv = MembraneConstants.PDBTM_FRAGMENTS_TM_BY_SEQUENCE_MOTIF_PATH.resolve("interactions-"
                + sequenceMotifDefinition.name() + ".tsv");

        long numberOfOccurrences = MembraneConstants.lines(interactionTsv)
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> split[0])
                .distinct()
                .count();
        int maximalInteractionIndex = Integer.valueOf(sequenceMotifDefinition.name().substring(2)) + 1;

        // initialize frequency table - for GG4 e.g. 4x4 matrix - non-symmetric (provides information on the direction)
        double[][] frequencyTable = new double[maximalInteractionIndex][];
        for(int index = 0; index < maximalInteractionIndex; index++) {
            frequencyTable[index] = new double[maximalInteractionIndex];
        }

        MembraneConstants.lines(interactionTsv)
                // filter header line
                .filter(line -> !line.startsWith("id"))
                // filter empty observations
                .filter(line -> !line.endsWith("-"))
                .map(line -> line.split("\t"))
                .map(split -> new Pair<>(Integer.valueOf(split[2]), Integer.valueOf(split[3])))
                .forEach(pair -> frequencyTable[pair.getLeft()][pair.getRight()] += (double) 1 / numberOfOccurrences);

        MembraneConstants.write(interactionTsv.getParent().resolve("interaction-frequency-" + sequenceMotifDefinition.name() + ".tsv"),
                Arrays.deepToString(frequencyTable).replace("],", "]," + System.lineSeparator()));
    }
}

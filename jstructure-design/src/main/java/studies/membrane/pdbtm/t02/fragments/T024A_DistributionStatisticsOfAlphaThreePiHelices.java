package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.Pair;
import studies.StudyConstants;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Compose one concise file of secondary structure information.
 *
 * Results:
 * Longer sequence motifs are more prone to disruptions (obviously). Especially motifs with glycine or proline have a
 * tendency to be part of unusual helix conformations. Maybe it is best to revert the approach and look for three-ten-
 * and pi-helices in the data set and find a pattern from that direction.
 *
 * Created by bittrich on 7/5/17.
 */
public class T024A_DistributionStatisticsOfAlphaThreePiHelices {
    public static void main(String[] args) {
        String output = Stream.of(SequenceMotifDefinition.values())
                .map(motif -> new Pair<>(motif, StudyConstants.lines(T024_LinkSequenceMotifsToSecondaryStructureElements.OUTPUT_PATH)
                        .filter(line -> !line.startsWith("id"))
                        .map(line -> line.split("\t"))
                        .filter(split -> split[1].equals(motif.name()))
                        .collect(Collectors.toList())))
                .map(T024A_DistributionStatisticsOfAlphaThreePiHelices::mapToLine)
                .collect(Collectors.joining(System.lineSeparator(),
                        "motif\tobservations\taminoAcids\tobsAlpha\tobsThree\tobsPi\tobsCoil\taaAlpha\taaThree\t" +
                                "aaPi\taaCoil\tobsImperfect\tobsBend\taaBend" + System.lineSeparator(),
                        ""));

        StudyConstants.write(T024_LinkSequenceMotifsToSecondaryStructureElements.OUTPUT_PATH.getParent().resolve("sse-by-sequence-motif.tsv"), output);
    }

    private static String mapToLine(Pair<SequenceMotifDefinition, List<String[]>> pair) {
        SequenceMotifDefinition motif = pair.getLeft();
        List<String[]> lines = pair.getRight();
        // number of observations
        int[] counts = new int[] { lines.size(),
                // number of amino acids
                lines.size() * Integer.valueOf(lines.get(0)[6]),
                // number of observations with alpha, three, pi or coil
                0, 0, 0, 0,
                // number of amino acids with alpha, three, pi or coil
                0, 0, 0, 0,
                // imperfect observations
                0,
                // number of observations with bend
                0,
                // number of amino acids with bend
                0 };
        for(String[] line : lines) {
            int alphaCount = Integer.valueOf(line[7]);
            int threeCount = Integer.valueOf(line[8]);
            int piCount = Integer.valueOf(line[9]);
            int coilCount = Integer.valueOf(line[10]);

            counts[2] += (alphaCount > 1) ? 1 : 0;
            counts[3] += (threeCount > 1) ? 1 : 0;
            counts[4] += (piCount > 1) ? 1 : 0;
            counts[5] += (coilCount > 1) ? 1 : 0;

            counts[6] += alphaCount;
            counts[7] += threeCount;
            counts[8] += piCount;
            counts[9] += coilCount;

            counts[10] += line[11].equals("true") ? 1 : 0;

            boolean hasBend = line[4].equals("true");
            counts[11] += hasBend ? 1 : 0;
            counts[12] += hasBend ? (line[5].split(", ").length + 1) : 0;
        }

        return motif.name() + "\t" + IntStream.of(counts)
                .mapToObj(Integer::toString)
                .collect(Collectors.joining("\t"));
    }
}

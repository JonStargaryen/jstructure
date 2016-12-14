package design.visualization;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.Pair;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Takes defined fragment clusters and moves them to the PyMOL directory, to visualize them later on.
 * Created by S on 13.12.2016.
 */
public class S01_MoveToPyMOLDirectory {
    private static final List<Pair<String, SequenceMotifDefinition>> configurations = Stream.of("tm-LY6", "trans-GG7", "trans-VL4")
            .map(line -> line.split("-"))
            .map(split -> new Pair<>(split[0], SequenceMotifDefinition.valueOf(split[1])))
            .collect(Collectors.toList());

    public static void main(String[] args) {
        configurations.forEach(System.out::println);
    }
}

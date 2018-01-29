package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.start2fold.Start2FoldConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class LowestQuartileSelector {
    private static final Logger logger = LoggerFactory.getLogger(LowestQuartileSelector.class);

    public LowestQuartileSelector(Path inputPath) {
        String[] filenameSplit = inputPath.toFile().getName().split("\\.");
        Path outputPath = inputPath.getParent().resolve(filenameSplit[0] + "-best." + filenameSplit[1]);

        List<String[]> lines;
        try(Stream<String> stream = Start2FoldConstants.lines(inputPath)) {
            lines = stream.map(line -> line.split(","))
                    .collect(Collectors.toList());
        }

        String[] header = lines.get(0);
        Map<String, List<String[]>> groupedLines = lines.stream()
                .skip(1)
                .collect(Collectors.groupingBy(split -> split[0] + "," + split[1] + "," + split[2]));
        // sort grouped lines by RMSD
        for(List<String[]> values : groupedLines.values()) {
            values.sort(Comparator.comparingDouble(value -> Double.valueOf(value[7])));
        }

        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        stringJoiner.add(Stream.of(header).collect(Collectors.joining(",")));
        for(List<String[]> values : groupedLines.values()) {
            int size = (int) (0.25 * values.size());
            logger.info("selecting {} of {} models in this bin",
                    size,
                    values.size());
            for(String[] value : values.subList(0, size)) {
                stringJoiner.add(Stream.of(value).collect(Collectors.joining(",")));
            }
        }

        Start2FoldConstants.write(outputPath, stringJoiner.toString());
    }
}

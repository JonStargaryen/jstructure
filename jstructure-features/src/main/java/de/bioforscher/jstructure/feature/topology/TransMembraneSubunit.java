package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.jstructure.model.structure.ParsingException;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * The generic instance of trans-membrane chain.
 * Created by bittrich on 5/17/17.
 */
public class TransMembraneSubunit {
    private String chainId;
    private String tilt;
    private List<IntegerRange> segments;
    private static final Pattern SEGMENTS_PATTERN = Pattern.compile(",");
    private static final Pattern RANGE_PATTERN = Pattern.compile("\\d+");

    TransMembraneSubunit(String describingText) {
        String[] globalSplit = describingText.split("-");
        this.chainId = globalSplit[0].trim();
        this.tilt = globalSplit[1].trim().replace("Tilt: ", "").trim();

        String segments = describingText.split("Segments:")[1].trim();
        this.segments = SEGMENTS_PATTERN.splitAsStream(segments)
                .map(this::mapToIntegerRange)
                .collect(Collectors.toList());
    }

    private IntegerRange mapToIntegerRange(String segment) {
        Matcher matcher = RANGE_PATTERN.matcher(segment);
        int start = Integer.MIN_VALUE;
        int end = Integer.MIN_VALUE;
        int counter = 0;
        while(matcher.find()) {
            if(counter == 1) {
                start = Integer.valueOf(matcher.group());
            }
            if(counter == 2) {
                end = Integer.valueOf(matcher.group());
            }
            counter++;
        }
        if(start == Integer.MIN_VALUE || end == Integer.MIN_VALUE) {
            throw new ParsingException("OPM output malformed for segment entry: " + segment);
        }
        return new IntegerRange(start, end);

    }

    public String getChainId() {
        return chainId;
    }

    public String getTilt() {
        return tilt;
    }

    public List<IntegerRange> getSegments() {
        return segments;
    }
}

package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.structure.selection.IntegerRange;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * The generic instance of trans-membrane helices.
 * Created by bittrich on 5/17/17.
 */
public class TransMembraneHelix {
    private String chainId;
    private String tilt;
    private List<IntegerRange> segments;
    private static final Pattern SEGMENTS_PATTERN = Pattern.compile(",");

    TransMembraneHelix(String describingText) {
        String[] globalSplit = describingText.split("-");
        this.chainId = globalSplit[0].trim();
        this.tilt = globalSplit[1].trim().replace("Tilt: ", "").trim();

        String segments = describingText.split("Segments:")[1].trim();
        this.segments = SEGMENTS_PATTERN.splitAsStream(segments)
                .map(rawString -> rawString.split("\\(")[1].split("\\)")[0])
                .map(rangeString -> rangeString.split("-"))
                .map(rangeSplit -> new IntegerRange(Integer.valueOf(rangeSplit[0].trim()), Integer.valueOf(rangeSplit[1].trim())))
                .collect(Collectors.toList());
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

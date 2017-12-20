package de.bioforscher.jstructure.contacts.visualization;

import java.util.List;

public class JsonModule {
    private final String identifier;
    private final List<JsonRange> ranges;

    public JsonModule(String identifier, List<JsonRange> ranges) {
        this.identifier = identifier;
        this.ranges = ranges;
    }

    public String getIdentifier() {
        return identifier;
    }

    public List<JsonRange> getRanges() {
        return ranges;
    }
}
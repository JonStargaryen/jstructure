package de.bioforscher.jstructure.contacts.visualization;

public class JsonSecondaryStructure extends JsonRange {
    private final String description;

    public JsonSecondaryStructure(int start, int end, String description) {
        super(start, end);
        this.description = description;
    }

    public String getDescription() {
        return description;
    }
}

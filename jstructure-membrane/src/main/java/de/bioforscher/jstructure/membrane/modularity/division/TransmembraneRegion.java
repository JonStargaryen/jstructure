package de.bioforscher.jstructure.membrane.modularity.division;

import org.jsoup.nodes.Element;

public class TransmembraneRegion {
    private final int start;
    private final int end;
    private final Type type;

    public enum Type {
        EXTRACELLULAR,
        HELICAL,
        CYTOPLASMIC
    }

    TransmembraneRegion(Element describingElement) {
        this.type = mapType(describingElement.attr("description"));
        this.start = Integer.valueOf(describingElement.getElementsByTag("begin").first().attr("position"));
        this.end = Integer.valueOf(describingElement.getElementsByTag("end").first().attr("position"));
    }

    private Type mapType(String description) {
        if(description.equals("Extracellular")) {
            return Type.EXTRACELLULAR;
        }
        if(description.equals("Cytoplasmic")) {
            return Type.CYTOPLASMIC;
        }
        if(description.startsWith("Helical")) {
            return Type.HELICAL;
        }

        throw new IllegalArgumentException("unexpected topology type: " + description);
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public Type getType() {
        return type;
    }

    @Override
    public String toString() {
        return start + " - " + end + " : " + type;
    }
}

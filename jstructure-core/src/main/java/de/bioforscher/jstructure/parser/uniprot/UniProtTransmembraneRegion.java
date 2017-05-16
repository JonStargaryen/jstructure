package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.parser.ParsingException;
import org.jsoup.nodes.Element;

/**
 * UniProt's topology annotation for membrane proteins.
 * Created by bittrich on 4/18/17.
 */
public class UniProtTransmembraneRegion {
    private int start, end;
    private Type type;

    public enum Type {
        EXTRACELLULAR,
        HELICAL,
        CYTOPLASMIC
    }

    public UniProtTransmembraneRegion() {
    }

    UniProtTransmembraneRegion(Element describingElement) {
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

        throw new ParsingException("unexpected topology type: " + description);
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
}

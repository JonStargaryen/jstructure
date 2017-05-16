package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.parser.ParsingException;
import org.jsoup.nodes.Element;

/**
 * UniProt's take on secondary structure elements.
 * Created by bittrich on 4/18/17.
 */
public class UniProtSecondaryStructureElement {
    private int start, end;
    private Type type;

    public enum Type {
        STRAND,
        HELIX
    }

    public UniProtSecondaryStructureElement() {
    }

    UniProtSecondaryStructureElement(Element describingElement) {
        this.start = Integer.valueOf(describingElement.getElementsByTag("begin").first().attr("position"));
        this.end = Integer.valueOf(describingElement.getElementsByTag("end").first().attr("position"));
        this.type = mapType(describingElement.attr("type"));
    }

    private Type mapType(String type) {
        if(type.equals("strand")) {
            return Type.STRAND;
        }
        if(type.equals("helix")) {
            return Type.HELIX;
        }

        throw new ParsingException("unexpected secondary structure element type: " + type);
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

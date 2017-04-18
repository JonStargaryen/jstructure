package de.bioforscher.jstructure.parser.uniprot;

import org.jsoup.nodes.Element;

/**
 * UniProt's annotation of disulfide bonds.
 * Created by bittrich on 4/18/17.
 */
public class UniProtDisulfideBond {
    private int start, end;

    public UniProtDisulfideBond() {
    }

    UniProtDisulfideBond(Element describingElement) {
        this.start = Integer.valueOf(describingElement.getElementsByTag("begin").first().attr("position"));
        this.end = Integer.valueOf(describingElement.getElementsByTag("end").first().attr("position"));
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }
}

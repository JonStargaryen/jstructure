package de.bioforscher.jstructure.feature.uniprot.old;

import org.jsoup.nodes.Element;

/**
 * UniProt's annotation of disulfide bonds.
 * Created by bittrich on 4/18/17.
 */
@Deprecated
public class UniProtDisulfideBond {
    private int start;
    private int end;

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

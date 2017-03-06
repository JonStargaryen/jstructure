package de.bioforscher.jstructure.parser.uniprot;

import org.jsoup.nodes.Element;

/**
 * Represents a UniProt citation.
 * Created by bittrich on 3/2/17.
 */
public class UniProtCitation {
    private String type, name, date, volume, first, last;

    UniProtCitation() {

    }

    UniProtCitation(Element describingElement) {
        this(describingElement.attr("type"),
                describingElement.attr("date"),
                describingElement.attr("name"),
                describingElement.attr("volume"),
                describingElement.attr("first"),
                describingElement.attr("last"));
    }

    public UniProtCitation(String type, String date, String name, String volume, String first, String last) {
        this.type = type;
        this.date = date;
        this.name = name;
        this.volume = volume;
        this.first = first;
        this.last = last;
    }

    public String getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    public String getDate() {
        return date;
    }

    public String getVolume() {
        return volume;
    }

    public String getFirst() {
        return first;
    }

    public String getLast() {
        return last;
    }
}
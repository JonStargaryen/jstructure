package de.bioforscher.jstructure.parser.uniprot;

import org.jsoup.nodes.Element;

/**
 * Annotation of active sites.
 * Created by bittrich on 4/6/17.
 */
public class UniProtActiveSite {
    private String position, description;

    UniProtActiveSite() {

    }

    UniProtActiveSite(Element describingElement) {
        this.position = describingElement.getElementsByAttribute("position").first().attr("position");
        this.description = describingElement.attr("description");
    }

    public String getPosition() {
        return position;
    }

    public String getDescription() {
        return description;
    }
}

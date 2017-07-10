package de.bioforscher.jstructure.feature.uniprot.old;

import org.jsoup.nodes.Element;

/**
 * Annotation of active sites.
 * Created by bittrich on 4/6/17.
 */
@Deprecated
public class UniProtActiveSite {
    private String position;
    private String description;

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

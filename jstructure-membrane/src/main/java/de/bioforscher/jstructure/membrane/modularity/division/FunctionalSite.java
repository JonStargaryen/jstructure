package de.bioforscher.jstructure.membrane.modularity.division;

import org.jsoup.nodes.Element;

public class FunctionalSite {
    private String position;
    private String description;

    FunctionalSite(Element describingElement) {
        this.position = describingElement.getElementsByAttribute("position").first().attr("position");
        this.description = describingElement.attr("description");
    }

    public String getPosition() {
        return position;
    }

    public String getDescription() {
        return description;
    }

    @Override
    public String toString() {
        return position + " : " + description;
    }
}

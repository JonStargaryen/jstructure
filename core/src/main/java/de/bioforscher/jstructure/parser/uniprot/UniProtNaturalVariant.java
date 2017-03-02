package de.bioforscher.jstructure.parser.uniprot;

import org.jsoup.nodes.Element;

/**
 * Describes a UniProt natural variant.
 * Created by bittrich on 3/2/17.
 */
public class UniProtNaturalVariant {
    private String id, original, variation, description, position;

    UniProtNaturalVariant() {

    }

    UniProtNaturalVariant(Element describingElement) {
        this(describingElement.attr("id"),
                describingElement.getElementsByTag("original").text(),
                describingElement.getElementsByTag("variation").text(),
                describingElement.attr("description"),
                describingElement.getElementsByTag("position").text());
    }

    public UniProtNaturalVariant(String id, String original, String variation, String description, String position) {
        this.id = id;
        this.original = original;
        this.variation = variation;
        this.description = description;
        this.position = position;
    }

    public String getId() {
        return id;
    }

    public String getOriginal() {
        return original;
    }

    public String getVariation() {
        return variation;
    }

    public String getDescription() {
        return description;
    }

    public String getPosition() {
        return position;
    }
}

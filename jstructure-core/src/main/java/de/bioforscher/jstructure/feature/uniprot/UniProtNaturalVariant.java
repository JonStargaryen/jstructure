package de.bioforscher.jstructure.feature.uniprot;

import org.jsoup.nodes.Element;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Describes a UniProt natural variant.
 * Created by bittrich on 3/2/17.
 */
public class UniProtNaturalVariant {
    private String id;
    private String original;
    private String variation;
    private String description;
    private List<String> evidence;
    private List<String> position;

    UniProtNaturalVariant() {

    }

    UniProtNaturalVariant(Element describingElement) {
        this(describingElement.attr("id"),
                describingElement.getElementsByTag("original").text(),
                describingElement.getElementsByTag("variation").text(),
                describingElement.attr("description"),
                describingElement.getElementsByTag("location").first().children().stream()
                        .map(element -> element.attr("position"))
                        .collect(Collectors.toList()),
                Pattern.compile("\\s").splitAsStream(describingElement.attr("evidence")).collect(Collectors.toList()));
    }

    public UniProtNaturalVariant(String id, String original, String variation, String description, List<String> position, List<String> evidence) {
        this.id = id;
        this.original = original;
        this.variation = variation;
        this.description = description;
        this.position = position;
        this.evidence = evidence;
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

    public List<String> getPosition() {
        return position;
    }

    public List<String> getEvidence() {
        return evidence;
    }
}

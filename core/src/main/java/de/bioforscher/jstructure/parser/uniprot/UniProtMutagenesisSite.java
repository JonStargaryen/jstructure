package de.bioforscher.jstructure.parser.uniprot;

import org.jsoup.nodes.Element;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Represents a UniProt mutation site and diseases linked to it.
 * Created by bittrich on 3/2/17.
 */
public class UniProtMutagenesisSite {
    private String original, variation, description;
    private List<String> evidence, position;

    UniProtMutagenesisSite() {

    }

    UniProtMutagenesisSite(Element describingElement) {
        this(describingElement.getElementsByTag("original").text(),
                describingElement.getElementsByTag("variation").text(),
                describingElement.attr("description"),
                describingElement.getElementsByTag("location").first().children().stream()
                        .map(element -> element.attr("position"))
                        .collect(Collectors.toList()),
                Pattern.compile("\\s").splitAsStream(describingElement.attr("evidence")).collect(Collectors.toList()));
    }

    public UniProtMutagenesisSite(String original, String variation, String description, List<String> position, List<String> evidence) {
        this.original = original;
        this.variation = variation;
        this.description = description;
        this.position = position;
        this.evidence = evidence;
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

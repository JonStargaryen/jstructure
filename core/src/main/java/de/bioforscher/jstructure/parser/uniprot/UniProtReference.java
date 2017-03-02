package de.bioforscher.jstructure.parser.uniprot;

import org.jsoup.nodes.Element;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a UniProt reference.
 * Created by bittrich on 3/2/17.
 */
public class UniProtReference {
    private String title, pubmed, doi;
    private List<String> authors;
    private UniProtCitation citation;

    UniProtReference() {

    }

    UniProtReference(Element describingElement) {
        this(describingElement.getElementsByTag("title").text(),
                describingElement.getElementsByAttributeValue("type", "PubMed").text(),
                describingElement.getElementsByAttributeValue("type", "DOI").text(),
                describingElement.getElementsByTag("person").stream()
                        .map(element -> element.attr("name"))
                        .collect(Collectors.toList()),
                new UniProtCitation(describingElement.getElementsByTag("citation").first()));
    }

    public UniProtReference(String title, String pubmed, String doi, List<String> authors, UniProtCitation citation) {
        this.title = title;
        this.pubmed = pubmed;
        this.doi = doi;
        this.authors = authors;
        this.citation = citation;
    }

    public String getTitle() {
        return title;
    }

    public String getPubmed() {
        return pubmed;
    }

    public String getDoi() {
        return doi;
    }

    public List<String> getAuthors() {
        return authors;
    }

    public UniProtCitation getCitation() {
        return citation;
    }
}

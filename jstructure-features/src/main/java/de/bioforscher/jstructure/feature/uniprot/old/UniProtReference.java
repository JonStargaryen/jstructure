package de.bioforscher.jstructure.feature.uniprot.old;

import org.jsoup.nodes.Element;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a UniProt reference.
 * Created by bittrich on 3/2/17.
 */
@Deprecated
public class UniProtReference {
    private String title;
    private String pubmed;
    private String doi;
    private List<String> authors;
    private UniProtCitation citation;

    UniProtReference() {

    }

    UniProtReference(Element describingElement) {
        try {
            this.title = describingElement.getElementsByTag("title").text();
        } catch (Exception e) {
            this.title = "";
        }
        try {
            this.pubmed = describingElement.getElementsByAttributeValue("type", "PubMed").attr("id");
        } catch (Exception e) {
            this.pubmed = "";
        }
        try {
            this.doi = describingElement.getElementsByAttributeValue("type", "DOI").attr("id");
        } catch (Exception e) {
            this.doi = "";
        }
        try {
            this.authors = describingElement.getElementsByTag("person").stream()
                    .map(element -> element.attr("name"))
                    .collect(Collectors.toList());
        } catch (Exception e) {
            this.authors = new ArrayList<>();
        }
        try {
            this.citation = new UniProtCitation(describingElement.getElementsByTag("citation").first());
        } catch (Exception e) {
            this.citation = null;
        }
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

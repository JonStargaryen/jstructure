package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.structure.Chain;
import org.jsoup.nodes.Document;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Gathers all UniProt information of a {@link Chain}.
 * Created by bittrich on 3/2/17.
 */
public class UniProtAnnotationContainer {
    private String uniProtId, sequence;
    private List<UniProtReference> references;
    private List<UniProtMutagenesisSite> mutagenesisSites;
    private List<UniProtNaturalVariant> naturalVariants;
    private List<UniProtActiveSite> activeSites;
    //TODO motifs?

    public UniProtAnnotationContainer() {
        this.uniProtId = "?";
        this.sequence = "?";
        this.references = new ArrayList<>();
        this.mutagenesisSites = new ArrayList<>();
        this.naturalVariants = new ArrayList<>();
        this.activeSites = new ArrayList<>();
    }

    UniProtAnnotationContainer(String uniProtId, Document describingDocument) {
        this.uniProtId = uniProtId;
        this.sequence = describingDocument.getElementsByTag("sequence").text().replaceAll("\\s+", "");
        this.references = describingDocument.getElementsByTag("reference").stream()
                .filter(element -> !element.getElementsByTag("citation").first().attr("type").equals("submission"))
                .map(UniProtReference::new)
                .collect(Collectors.toList());
        this.mutagenesisSites = describingDocument.getElementsByAttributeValue("type", "mutagenesis site").stream()
                .map(UniProtMutagenesisSite::new)
                .collect(Collectors.toList());
        this.naturalVariants = describingDocument.getElementsByAttributeValue("type", "sequence variant").stream()
                .map(UniProtNaturalVariant::new)
                .collect(Collectors.toList());
        this.activeSites = describingDocument.getElementsByAttributeValue("type", "active site").stream()
                .map(UniProtActiveSite::new)
                .collect(Collectors.toList());
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public String getSequence() {
        return sequence;
    }

    public List<UniProtReference> getReferences() {
        return references;
    }

    public List<UniProtMutagenesisSite> getMutagenesisSites() {
        return mutagenesisSites;
    }

    public List<UniProtNaturalVariant> getNaturalVariants() {
        return naturalVariants;
    }

    public List<UniProtActiveSite> getActiveSites() {
        return activeSites;
    }
}

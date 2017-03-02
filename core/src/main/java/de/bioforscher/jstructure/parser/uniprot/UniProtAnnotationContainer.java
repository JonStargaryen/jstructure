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
    private List<UniProtReference> references;
    private List<UniProtMutagenesisSite> mutagenesisSites;
    private List<UniProtNaturalVariant> naturalVariants;
    //TODO motifs?

    public UniProtAnnotationContainer() {
        this.references = new ArrayList<>();
        this.mutagenesisSites = new ArrayList<>();
        this.naturalVariants = new ArrayList<>();
    }

    UniProtAnnotationContainer(Document describingDocument) {
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
}

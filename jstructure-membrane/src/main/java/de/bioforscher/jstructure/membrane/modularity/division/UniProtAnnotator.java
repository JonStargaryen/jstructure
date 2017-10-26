package de.bioforscher.jstructure.membrane.modularity.division;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;

import java.io.IOException;
import java.io.UncheckedIOException;

public class UniProtAnnotator {
    private static final String UNIPROT_FETCH_URL = "http://www.uniprot.org/uniprot/%s.xml";

    public UniProtAnnotationContainer process(String uniProtId) {
        try {
            Document document = Jsoup.connect(String.format(UNIPROT_FETCH_URL, uniProtId)).get();
            return new UniProtAnnotationContainer(document);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}

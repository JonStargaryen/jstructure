package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;

import static de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator.UNIPROT_ANNOTATION;

/**
 * Queries the UniProt db and annotates available information to each protein chain.
 * Created by bittrich on 3/2/17.
 */
@FeatureProvider(provides = { UNIPROT_ANNOTATION }, requires = {SiftsParser.UNIPROT_ID })
public class UniProtAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(UniProtAnnotator.class);
    public static final String UNIPROT_ANNOTATION = "UNIPROT_ANNOTATION";
    private static final String UNIPROT_FETCH_URL = "http://www.uniprot.org/uniprot/%s.xml";

    @Override
    protected void processInternally(Protein protein) {
        protein.chains()
                // for safety: ignore non-amino-acid chains
                .filter(chain -> chain.aminoAcids().count() > 0)
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        String uniprotId = chain.getFeature(String.class, SiftsParser.UNIPROT_ID);

        if(uniprotId.equals(SiftsParser.UNKNOWN_MAPPING)) {
            chain.setFeature(UNIPROT_ANNOTATION, new UniProtAnnotationContainer());
            return;
        }

        try {
            logger.debug("fetching UniProt information for {}", uniprotId);
            //TODO error-handling for connectivity issues and downtime
            UniProtAnnotationContainer uniProtAnnotationContainer = new UniProtAnnotationContainer(Jsoup.connect(String.format(UNIPROT_FETCH_URL, uniprotId)).get());
            chain.setFeature(UNIPROT_ANNOTATION, uniProtAnnotationContainer);
            logger.info("{} mutations, {} variants, {} references", uniProtAnnotationContainer.getMutagenesisSites().size(), uniProtAnnotationContainer.getNaturalVariants().size(), uniProtAnnotationContainer.getReferences().size());
        } catch (IOException e) {
            throw new UncheckedIOException("could not retrieve UniProt entry " + uniprotId + ": " + e.getCause().getLocalizedMessage(), e);
        }
    }
}

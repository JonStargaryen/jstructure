package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.sifts.ChainSiftsMapping;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;
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
@FeatureProvider(provides = { UNIPROT_ANNOTATION }, requires = { SiftsMappingProvider.SIFTS_MAPPING })
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
        String uniprotId = chain.getFeature(ChainSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING).getUniProtId();

        UniProtAnnotationContainer uniProtAnnotationContainer = process(uniprotId);
        chain.setFeature(UNIPROT_ANNOTATION, uniProtAnnotationContainer);
        logger.info("{} mutations, {} variants, {} references", uniProtAnnotationContainer.getMutagenesisSites().size(), uniProtAnnotationContainer.getNaturalVariants().size(), uniProtAnnotationContainer.getReferences().size());
    }

    public UniProtAnnotationContainer process(String uniprotId) {
        //TODO error-handling for connectivity issues and downtime
        try {
            if(uniprotId.equals(SiftsMappingProvider.UNKNOWN_MAPPING)) {
                return new UniProtAnnotationContainer();
            }

            logger.debug("fetching UniProt information for {}", uniprotId);
            return new UniProtAnnotationContainer(uniprotId, Jsoup.connect(String.format(UNIPROT_FETCH_URL, uniprotId)).get());
        } catch (IOException e) {
            throw new UncheckedIOException("could not retrieve UniProt entry " + uniprotId + ": " + e.getCause().getLocalizedMessage(), e);
        }
    }
}

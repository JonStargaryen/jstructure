package de.bioforscher.jstructure.feature.uniprot.old;

import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.feature.mapping.ChainMapping;
import de.bioforscher.jstructure.feature.mapping.ResidueMapping;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * Queries the UniProt db and annotates available information to each protein chain.
 * Created by bittrich on 3/2/17.
 */
@FeatureProvider(provides = UniProtAnnotationContainer.class, requires = ResidueMapping.class)
@Deprecated
public class UniProtAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(UniProtAnnotator.class);
    private static final String UNIPROT_FETCH_URL = "http://www.uniprot.org/uniprot/%s.xml";

    @Override
    protected void processInternally(Protein protein) {
        protein.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        String uniprotId = chain.getFeatureContainer().getFeature(ChainMapping.class).getUniProtId();
        // mapping may fail
        if(ChainMapping.UNKNOWN_MAPPING.equals(uniprotId)) {
            logger.debug("could not retrieve UniProt mapping for {} - ignoring chain", chain.getChainIdentifier().getFullName());
            return;
        }

        UniProtAnnotationContainer uniProtAnnotationContainer = processUniProtId(uniprotId);
        chain.getFeatureContainer().addFeature(uniProtAnnotationContainer);
        logger.debug("{} mutations, {} variants, {} references", uniProtAnnotationContainer.getMutagenesisSites().size(), uniProtAnnotationContainer.getNaturalVariants().size(), uniProtAnnotationContainer.getReferences().size());
    }

    public UniProtAnnotationContainer processUniProtId(String uniProtId) {
        try {
            logger.debug("fetching UniProt information for {}", uniProtId);
            return new UniProtAnnotationContainer(this, uniProtId,
                    Jsoup.connect(String.format(UNIPROT_FETCH_URL, uniProtId)).get());
        } catch (IOException e) {
            throw new ComputationException("could not retrieve UniProt entry " + uniProtId + ": " +
                    e.getCause().getLocalizedMessage(), e);
        }
    }
}

package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.NoSuchElementException;
import java.util.Optional;


/**
 * Maps PDB identifiers to UniProt.
 * Created by bittrich on 3/1/17.
 */

@FeatureProvider(provides = { SiftsParser.UNIPROT_ID }, priority = 50)
public class UniProtMapper extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(UniProtMapper.class);
    private static final String MAPPING_REST_URL = "http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&&id=%s&chain=%s";

    @Override
    protected void processInternally(Protein protein) {
        protein.chains().forEach(this::processInternally);
    }

    protected void processInternally(Chain chain) {
        String pdbId = chain.getParentProtein().getName();
        try {
            String uniProtId = parseUniprotInformation(new URL(String.format(MAPPING_REST_URL, pdbId, chain.getChainId())).openStream()).orElseThrow(NoSuchElementException::new);
            logger.debug("mapped UniProt id for {} is {}", pdbId, uniProtId);
            chain.setFeature(SiftsParser.UNIPROT_ID, uniProtId);
        } catch (IOException | NoSuchElementException e) {
            logger.warn("could not map to UniProt for {}-{}", chain.getParentProtein().getName(), chain.getChainId());
            chain.setFeature(SiftsParser.UNIPROT_ID, SiftsParser.UNKNOWN_MAPPING);
        }
    }

    private Optional<String> parseUniprotInformation(InputStream inputStream) throws IOException {
        try (InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
            try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                return bufferedReader.lines()
                        .filter(line -> line.startsWith("AC"))
                        .map(line -> line.split("AC: ")[1])
                        .findFirst();
            }
        }
    }
}

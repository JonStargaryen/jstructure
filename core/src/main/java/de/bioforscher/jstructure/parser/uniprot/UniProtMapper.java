package de.bioforscher.jstructure.parser.uniprot;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ParsingException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.NoSuchElementException;
import java.util.Optional;

import static de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator.UNIPROT_ID;

/**
 * Maps PDB identifiers to UniProt.
 * Created by bittrich on 3/1/17.
 */
class UniProtMapper {
    private static final Logger logger = LoggerFactory.getLogger(UniProtMapper.class);
    private static final String MAPPING_REST_URL = "http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&&id=%s&chain=%s";

    public static void process(Protein protein) {
        protein.chains().forEach(UniProtMapper::process);
    }

    public static void process(Chain chain) {
        String pdbId = chain.getParentProtein().getName();
        try {
            String uniProtId = parseUniprotInformation(new URL(String.format(MAPPING_REST_URL, pdbId, chain.getChainId())).openStream()).orElseThrow(NoSuchElementException::new);
            logger.debug("mapped UniProt id for {} is {}", pdbId, uniProtId);
            chain.setFeature(UNIPROT_ID, uniProtId);
        } catch (IOException | NoSuchElementException e) {
            throw new ParsingException("could not determine UniProt mapping for " + pdbId + " - " + chain.getChainId(), e);
        }
    }

    private static Optional<String> parseUniprotInformation(InputStream inputStream) throws IOException {
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

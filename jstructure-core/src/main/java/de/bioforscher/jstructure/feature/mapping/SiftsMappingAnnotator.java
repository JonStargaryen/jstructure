package de.bioforscher.jstructure.feature.mapping;

import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.parser.Parser;
import org.jsoup.select.Elements;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.net.URL;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

/**
 * Retrieves cross-database mappings for residues and chains using the Sifts-project.
 * Created by bittrich on 5/17/17.
 */
@FeatureProvider(provides = { ResidueMapping.class, ChainMapping.class })
public class SiftsMappingAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(SiftsMappingAnnotator.class);

    // the SIFTS-web-resource
    private static final String XML_FETCH_URL = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/%s/%s.xml.gz";
    // local resources
    private static final String SIFTS_DIR = "sifts/";
    private static final String ENZYME_MAPPING = SIFTS_DIR + "pdb_chain_enzyme.csv";
    private static final String PFAM_MAPPING = SIFTS_DIR + "pdb_chain_pfam.csv";

    private static final String UNKNOWN_MAPPING = "?";
    private static final String[] UNKNOWN_PFAM_MAPPING = new String[] { UNKNOWN_MAPPING, UNKNOWN_MAPPING, UNKNOWN_MAPPING, UNKNOWN_MAPPING };

    @Override
    protected void processInternally(Protein protein) {
        Document document = downloadXml(protein.getPdbId().getPdbId());
        protein.aminoAcidChains()
                .forEach(chain -> processInternally(document, chain));
    }

    private void processInternally(Document document, Chain chain) {
        chain.aminoAcids().forEach(group -> {
            ResidueMapping mapping = mapGroup(document, chain.getChainId().getChainId(), group.getResidueNumber());
            group.getFeatureContainer().addFeature(mapping);
        });

        String pdbId = chain.getParentProtein().getPdbId().getPdbId();
        String ecNumber = getLinesForPdbId(ENZYME_MAPPING, pdbId)
                .filter(split -> split[1].equals(chain.getChainId()))
                .map(split -> split[3])
                .findFirst()
                .orElse(UNKNOWN_MAPPING);

        String[] mappingString = getLinesForPdbId(PFAM_MAPPING, pdbId)
                .filter(split -> split[1].equals(chain.getChainId()))
                .findFirst()
                .orElse(UNKNOWN_PFAM_MAPPING);

        chain.getFeatureContainer().addFeature(new ChainMapping(this, mappingString[2], ecNumber, mappingString[3]));
    }

    private Stream<String> getLines(String path) {
        return getLinesFromResource(path)
                .filter(line -> !line.startsWith("#") && !line.startsWith("PDB"));
    }

    private Stream<String[]> getLinesForPdbId(String path, String pdbId) {
        return getLines(path)
                .filter(line -> line.startsWith(pdbId.toLowerCase()))
                .map(line -> line.split(","));
    }

    /**
     * Resolves the URL to the SIFTS-mapping file for a given pdb-id, fetches is, decompresses it (as it is delivered as
     * <code>gz</code>, parses its xml and returns its Jsoup-representation.
     * @param pdbId the id to fetch
     * @return the Jsoup-document describing the return XML-file
     */
    Document downloadXml(String pdbId) {
        if(pdbId.length() != 4) {
            throw new IllegalArgumentException(pdbId + " is no valid pdb-id");
        }

        try {
            pdbId = pdbId.toLowerCase();
            String url = String.format(XML_FETCH_URL, pdbId.substring(1, 3), pdbId);
            InputStream inputStream = new BufferedInputStream(new URL(url).openStream(), 65535);
            try(GZIPInputStream gzipInputStream = new GZIPInputStream(inputStream)) {
                return Jsoup.parse(gzipInputStream, "UTF-8", url, Parser.xmlParser());
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    ResidueMapping mapGroup(final Document document, final String chainId, final int pdbResidueNumber) {
        return mapGroup(document, chainId, String.valueOf(pdbResidueNumber));
    }

    private ResidueMapping mapGroup(final Document document, final String chainId, final String pdbResidueNumber) {
        Element describingElement = mapToDescribingElement(document, chainId, pdbResidueNumber);
        Elements uniProtElements = describingElement.getElementsByAttributeValue("dbSource", "UniProt");
        if(!uniProtElements.isEmpty()) {
            Element uniProtElement = uniProtElements.first();
            return new ResidueMapping(this, uniProtElement.attr("dbResNum"));
        } else {
            logger.warn("could not retrieve UniProt mapping for " + chainId + "-" + pdbResidueNumber);
            return new ResidueMapping(this);
        }
    }

    Element mapToDescribingElement(final Document document, final String chainId, final String pdbResidueNumber) {
        return document.getElementsByTag("crossRefDb").stream()
                // filter for PDB-related entries
                .filter(element -> element.attr("dbSource").equals("PDB"))
                // filter for residue number
                .filter(element -> element.attr("dbResNum").equals(pdbResidueNumber))
                // and chain
                .filter(element -> element.attr("dbChainId").equals(chainId))
                .map(Element::parent)
                .findFirst()
                .orElseThrow(() -> new ComputationException("could not obtain SIFTS-mapping for residue " + chainId + "-" + pdbResidueNumber));
    }
}

package aars;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Gather all SIFTS information for a given protein chain.
 * @link http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
 * Created by bittrich on 3/13/17.
 */
class SiftsParser extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(SiftsParser.class);
    private static final String SIFTS_DIR = "sifts/";
    public static final String UNIPROT_ID = "UNIPROT_ID";
    public static final String PFAM_ID = "PFAM_ID";
    public static final String EC_NUMBER = "EC_NUMBER";
    public static final String UNKNOWN_MAPPING = "?";
    private static final String[] UNKNOWN_PFAM_MAPPING = new String[] { UNKNOWN_MAPPING, UNKNOWN_MAPPING, UNKNOWN_MAPPING, UNKNOWN_MAPPING };
    private static final String ENZYME_MAPPING = SIFTS_DIR + "pdb_chain_enzyme.csv";
    private static final String PFAM_MAPPING = SIFTS_DIR + "pdb_chain_pfam.csv";

    @Override
    protected void processInternally(Protein protein) {
        protein.chains()
                // for safety: ignore non-amino-acid chains
                .filter(chain -> chain.aminoAcids().count() > 0)
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        String pdbId = chain.getParentProtein().getName().toLowerCase();
        chain.setFeature(EC_NUMBER, getLines(ENZYME_MAPPING, pdbId)
                .filter(split -> split[1].equals(chain.getChainId()))
                .map(split -> split[3])
                .findFirst()
                .orElse(UNKNOWN_MAPPING));

        String[] mappingString = getLines(PFAM_MAPPING, pdbId)
                .filter(split -> split[1].equals(chain.getChainId()))
                .findFirst()
                .orElse(UNKNOWN_PFAM_MAPPING);

        chain.setFeature(UNIPROT_ID, mappingString[2]);
        chain.setFeature(PFAM_ID, mappingString[3]);

        logger.debug("mapping of '{}': {}", chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId(), Arrays.toString(mappingString));
    }

    public List<String> mapToUniProt(String pdbId, String chainId) {
        return getLines(PFAM_MAPPING, pdbId)
                .filter(split -> split[1].equals(chainId))
                .map(split -> split[2])
                .distinct()
                .collect(Collectors.toList());
    }

    private Stream<String[]> getLines(String path, String pdbId) {
        return getLinesFromResource(path)
                .filter(line -> !line.startsWith("#") && !line.startsWith("PDB"))
                .filter(line -> line.startsWith(pdbId.toLowerCase()))
                .map(line -> line.split(","));
    }
}

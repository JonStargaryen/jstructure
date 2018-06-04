package de.bioforscher.jstructure.model.identifier;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.regex.Pattern;

/**
 * Convenience class to build identifiers.
 * Created by bittrich on 7/11/17.
 */
public class IdentifierFactory {
    private static final Logger logger = LoggerFactory.getLogger(IdentifierFactory.class);
    private static final Pattern EC_PATTERN = Pattern.compile("^([1-6])$|^([1-6])\\.(\\d{1,2})$|^([1-6])\\.(\\d{1,2})\\.(\\d{1,2})$|^([1-6])\\.(\\d{1,2})\\.(\\d{1,2})\\.(\\d{1,3})$");
    private static final Pattern PDBID_PATTERN = Pattern.compile("[1-9][A-Za-z0-9]{3}");
    private static final Pattern PFAM_PATTERN = Pattern.compile("PF\\d{5}");
    private static final Pattern UNIPROT_MODERN_PATTERN = Pattern.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}");

    public static EnzymeClassificationIdentifier createEnzymeClassificationIdentifier(String ecNumber) {
        validate(ecNumber, EC_PATTERN, "ecNumber");
        return new EnzymeClassificationIdentifier(ecNumber);
    }

    public static ProteinFamilyIdentifier createProteinFamilyIdentifier(String pfamId) {
        validate(pfamId, PFAM_PATTERN, "pfamId");
        return new ProteinFamilyIdentifier(pfamId);
    }

    public static UniProtIdentifier createUniProtIdentifier(String uniProtId) {
        // uniProtId matches new format, straight forward create the instance
        if(UNIPROT_MODERN_PATTERN.matcher(uniProtId).find()) {
            return new UniProtIdentifier(uniProtId);
        }

        // the hard case, let the id be resolved by UniProt
        try {
            InputStream inputStream = new URL("http://www.uniprot.org/uniprot/" + uniProtId + ".fasta").openStream();
            try (InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    // >sp|P38554|CYC32_DESNO Cytochrome c3, 26 kDa OS=Desulfomicrobium norvegicum (strain DSM 1741 / NCIMB 8310) PE=1 SV=1
                    String header = bufferedReader.readLine();
                    return new UniProtIdentifier(header.split("\\|")[1]);
                }
            }
        } catch (IOException e) {
            throw new IllegalArgumentException("'" + uniProtId + "' is no valid uniProtId", e);
        }
    }

    /**
     * This protein is identified by only a pdbId.
     * @see #createProteinIdentifier(String, String)
     */
    public static ProteinIdentifier createProteinIdentifier(String pdbId) {
        if(pdbId == null || pdbId.isEmpty()) {
            throw new IllegalArgumentException("pdbId cannot be empty or null, when no additional name is provided");
        }
        return createProteinIdentifier(pdbId, "");
    }

    /**
     * Construct an identifier of a pdbId and an additional name.
     * @param pdbId the pdbId
     * @param additionalName further information
     * @return the create instance
     * @throws IllegalArgumentException when pdbId does not match the pdbId pattern
     */
    public static ProteinIdentifier createProteinIdentifier(String pdbId, String additionalName) {
        if((pdbId == null || pdbId.isEmpty()) && (additionalName == null || additionalName.isEmpty())) {
            throw new IllegalArgumentException("pdbId and additionalName cannot be both empty or null");
        }
        if(pdbId != null && !pdbId.isEmpty()) {
            validate(pdbId, PDBID_PATTERN, "pdbId");
        }
        validateAdditionalName(additionalName);
        return new ProteinIdentifier(pdbId, additionalName);
    }

    /**
     * Creates a ChainIdentifier by a given pdbId and chainId.
     * @param pdbId the pdbId of the structure
     * @param chainId the chainId of this chain
     * @return the identifier
     */
    public static ChainIdentifier createChainIdentifier(String pdbId, String chainId) {
        return createChainIdentifier(createProteinIdentifier(pdbId), chainId);
    }

    /**
     * Creates a new instance of a ChainIdentifer.
     * @param proteinIdentifier the parent identifier
     * @param chainId this chain's id
     * @return the create instance
     */
    public static ChainIdentifier createChainIdentifier(ProteinIdentifier proteinIdentifier, String chainId) {
        return new ChainIdentifier(proteinIdentifier, chainId);
    }

    /**
     * Create a ResidueIdentifier merely by the residue number.
     * @param residueNumber the group's number
     * @return the identifier instance
     */
    public static ResidueIdentifier createResidueIdentifier(int residueNumber) {
        return createResidueIdentifier(residueNumber, "");
    }

    /**
     * Create a ResidueIdentifier with an addtional insertion code.
     * @param residueNumber the group's number
     * @param insertionCode the group's insertion code
     * @return the identifier instance
     */
    public static ResidueIdentifier createResidueIdentifier(int residueNumber, String insertionCode) {
        return new ResidueIdentifier(residueNumber, insertionCode);
    }

    private static void validate(String id, Pattern pattern, String name) {
        if(id == null || !pattern.matcher(id).find()) {
            throw new IllegalArgumentException("'" + id + "' is no valid " + name);
        }
    }

    private static void validateAdditionalName(String additionalName) {
        if(additionalName.contains("_")) {
//            throw new IllegalArgumentException("'" + additionalName + "' is no valid additional protein name as it " +
//                    "contains a underscore _ which is used to separate the chain identifer from the pdbId");
            logger.debug("'" + additionalName + "' is no valid additional protein name as it " +
                    "contains a underscore _ which is used to separate the chain identifer from the pdbId");
        }
    }

    /**
     * Parses residue identifiers with potential insertion codes. E.g., will handle '95' as well as '95A' and map them
     * to residue number 95 respectively reside number 95 with insertion code A.
     * @return the {@link ResidueIdentifier} instance
     */
    public static ResidueIdentifier createResidueIdentifier(String residueNumberString) {
        try {
            // plain number
            return createResidueIdentifier(Integer.valueOf(residueNumberString));
        } catch (NumberFormatException e) {
            // has insertion code
            int splitPoint = residueNumberString.length() - 1;
            return createResidueIdentifier(Integer.valueOf(residueNumberString.substring(0, splitPoint)),
                    residueNumberString.substring(splitPoint));
        }
    }
}

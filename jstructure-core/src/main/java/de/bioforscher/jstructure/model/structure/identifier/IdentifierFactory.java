package de.bioforscher.jstructure.model.structure.identifier;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.regex.Pattern;

/**
 * Convenience class to build identifiers.
 * Created by bittrich on 7/11/17.
 */
public class IdentifierFactory {
    private static final Logger logger = LoggerFactory.getLogger(IdentifierFactory.class);
    private static final Pattern PDBID_PATTERN = Pattern.compile("[1-9][A-Za-z0-9]{3}");

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
            validatePdbId(pdbId);
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


    private static void validatePdbId(String pdbId) {
        if(!PDBID_PATTERN.matcher(pdbId).find()) {
            throw new IllegalArgumentException("'" + pdbId + "' is no valid pdbId");
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
}

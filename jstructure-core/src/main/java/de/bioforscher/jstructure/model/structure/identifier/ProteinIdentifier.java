package de.bioforscher.jstructure.model.structure.identifier;

import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents the name/id of a protein. Either a standard PDB-id or a more complex name such as a filename or added
 * information. PDB-ids will internally always be represented by lowercase strings.
 * Created by bittrich on 4/27/17.
 */
public class ProteinIdentifier {
    private static final Pattern PDBID_PATTERN = Pattern.compile("[1-9][A-Za-z0-9]{3}");
    private final String pdbId;
    private final String additionalName;
    //TODO maybe some selection criteria and stuff automatically
    //TODO maybe the source / filename automatically
    //TODO consistent naming pattern

    public static final ProteinIdentifier UNKNOWN_PROTEIN_ID = ProteinIdentifier.createFromAdditionalName("UNKNOWN-PROTEIN");

    private ProteinIdentifier(String pdbId, String additionalName) {
        this.pdbId = pdbId.toLowerCase();
        this.additionalName = additionalName;
    }

    private static void validatePdbId(String pdbId) {
        if(!PDBID_PATTERN.matcher(pdbId).find()) {
            throw new IllegalArgumentException("'" + pdbId + "' is no valid pdbId");
        }
    }

    private static void validateAdditionalName(String additionalName) {
        if(additionalName.contains("_")) {
            throw new IllegalArgumentException("'" + additionalName + "' is no valid additional protein name as it " +
                    "contains a underscore _ which is used to separate the chain identifer from the pdbId");
        }
    }

    /**
     * This protein is identified by only a pdbId.
     * @see #createFromPdbIdAndName(String, String)
     */
    public static ProteinIdentifier createFromPdbId(String pdbId) {
        return createFromPdbIdAndName(pdbId, "");
    }

    /**
     * This protein is not identified by any pdbId, but a custom string instead.
     * @see #createFromPdbIdAndName(String, String)
     */
    public static ProteinIdentifier createFromAdditionalName(String additionalName) {
        return createFromPdbIdAndName("", additionalName);
    }

    /**
     * Construct an identifier of a pdbId and an additional name.
     * @param pdbId the pdbId
     * @param additionalName further information
     * @return the create instance
     * @throws IllegalArgumentException when pdbId does not match the pdbId pattern
     */
    public static ProteinIdentifier createFromPdbIdAndName(String pdbId, String additionalName) {
        if(!pdbId.isEmpty()) {
            validatePdbId(pdbId);
        }
        validateAdditionalName(additionalName);
        return new ProteinIdentifier(pdbId, additionalName);
    }

    /**
     * Returns the pdbId of this protein (if any).
     * @return the protein's pdbId
     */
    public String getPdbId() {
        return pdbId;
    }

    /**
     * Returns the additional naming string of this protein (if any).
     * @return the protein's additional name
     */
    public String getAdditionalName() {
        return additionalName;
    }

    /**
     * Composes the complex name of this protein (i.e. its pdbId and additional name, delimited by a hyphen).
     * @return the complex name of this protein
     */
    public String getFullName() {
        return Stream.of(pdbId, additionalName)
                .filter(string -> !string.isEmpty())
                .collect(Collectors.joining("-"));
    }

    /**
     * Checks whether this identifier refers to a standard pdbId and nothing else.
     * @return <code>true</code> iff this protein directly originates from a PDB-entry
     */
    public boolean isStandardPdbId() {
        // identifier can only by standard, when the pdbId is set and no additional name was provided
        return !pdbId.isEmpty() && additionalName.isEmpty();
    }

    @Override
    public String toString() {
        return getFullName();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ProteinIdentifier pdbId1 = (ProteinIdentifier) o;

        return (pdbId != null ? pdbId.equals(pdbId1.pdbId) : pdbId1.pdbId == null) && (additionalName != null ? additionalName.equals(pdbId1.additionalName) : pdbId1.additionalName == null);
    }

    @Override
    public int hashCode() {
        int result = pdbId != null ? pdbId.hashCode() : 0;
        result = 31 * result + (additionalName != null ? additionalName.hashCode() : 0);
        return result;
    }
}

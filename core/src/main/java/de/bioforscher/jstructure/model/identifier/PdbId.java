package de.bioforscher.jstructure.model.identifier;

import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents the name/id of a protein. Either a standard PDB-id or a more complex name such as a filename or added
 * information. PDB-ids will internally always be represented by lowercase strings.
 * Created by bittrich on 4/27/17.
 */
public class PdbId {
    private static final Pattern PDBID_PATTERN = Pattern.compile("[1-9][A-Za-z0-9]{3}");
    private String pdbId, additionalName;
    //TODO maybe some selection criteria and stuff automatically
    //TODO maybe the source / filename automatically

    private PdbId(String pdbId, String additionalName) {
        this.pdbId = pdbId.toLowerCase();
        this.additionalName = additionalName;
    }

    private static void validatePdbId(String pdbId) {
        if(!PDBID_PATTERN.matcher(pdbId).find()) {
            throw new IllegalArgumentException("'" + pdbId + "' is no valid PDB-id");
        }
    }

    /**
     * @see #createFromPdbIdAndName(String, String)
     */
    public static PdbId createFromPdbId(String pdbId) {
        validatePdbId(pdbId);
        return new PdbId(pdbId, "");
    }

    /**
     * @see #createFromPdbIdAndName(String, String)
     */
    public static PdbId createFromName(String additionalName) {
        return new PdbId("", additionalName);
    }

    /**
     * Construct an identifier of a PDB-id and an additional name.
     * @param pdbId the PDB-id
     * @param additionalName further information
     * @return the create instance
     * @throws IllegalArgumentException when pdbId does not match the PDB-id pattern
     */
    public static PdbId createFromPdbIdAndName(String pdbId, String additionalName) {
        validatePdbId(pdbId);
        return new PdbId(pdbId, additionalName);
    }

    /**
     * Returns the PDB-id of this protein (if any).
     * @return the protein's PDB-id
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
     * Composes the complex name of this protein (i.e. its PDB-id and additional name).
     * @return the complex name of this protein
     */
    public String getFullName() {
        return Stream.of(pdbId, additionalName)
                .filter(string -> !string.isEmpty())
                .collect(Collectors.joining("-"));
    }

    /**
     * Checks whether this identifier refers to a standard PDB-id and nothing else.
     * @return <code>true</code> iff this protein directly originates from a PDB-entry
     */
    public boolean isStandardPdbId() {
        return PDBID_PATTERN.matcher(pdbId).find() && additionalName.isEmpty();
    }

    @Override
    public String toString() {
        return getFullName();
    }
}

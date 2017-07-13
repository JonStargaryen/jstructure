package de.bioforscher.jstructure.model.structure.identifier;

import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents the name/id of a protein. Either a standard PDB-id or a more complex name such as a filename or added
 * information. PDB-ids will internally always be represented by lowercase strings. Instances of this class are
 * immutable.
 * Created by bittrich on 4/27/17.
 */
public class ProteinIdentifier {
    private final String pdbId;
    private final String additionalName;
    //TODO maybe some selection criteria and stuff automatically
    //TODO maybe the source / filename automatically
    //TODO consistent naming pattern

    public static final ProteinIdentifier UNKNOWN_PROTEIN_ID = IdentifierFactory.createProteinIdentifier("", "UNKNOWN-PROTEIN");

    ProteinIdentifier(String pdbId, String additionalName) {
        this.pdbId = pdbId.toLowerCase();
        this.additionalName = additionalName;
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
     * Checks whether this identifier refers to a standard pdbId and nothing else. Identifiers can only by standard,
     * when the pdbId is set and no additional name was provided.
     * @return <code>true</code> iff this protein directly originates from a PDB-entry
     */
    public boolean isStandardPdbId() {
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

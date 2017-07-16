package de.bioforscher.jstructure.model.identifier;

/**
 * Represents UniProt identifiers.
 * Created by S on 14.07.2017.
 */
public class UniProtIdentifier extends AbstractIdentifier {
    private final String uniProtId;

    UniProtIdentifier(String uniProtId) {
        this.uniProtId = uniProtId;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    @Override
    public String toString() {
        return uniProtId;
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;

        UniProtIdentifier that = (UniProtIdentifier) other;

        return uniProtId.equals(that.uniProtId);
    }

    @Override
    public int hashCode() {
        return uniProtId.hashCode();
    }
}

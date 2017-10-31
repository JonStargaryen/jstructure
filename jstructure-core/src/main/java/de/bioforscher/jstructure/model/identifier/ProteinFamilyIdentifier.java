package de.bioforscher.jstructure.model.identifier;

/**
 * The PFAM identifier.
 * Created by S on 14.07.2017.
 */
public class ProteinFamilyIdentifier extends AbstractIdentifier {
    private final String pfamId;

    ProteinFamilyIdentifier(String pfamId) {
        this.pfamId = pfamId;
    }

    public String getPfamId() {
        return pfamId;
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;

        ProteinFamilyIdentifier that = (ProteinFamilyIdentifier) other;

        return pfamId.equals(that.pfamId);
    }

    @Override
    public int hashCode() {
        return pfamId.hashCode();
    }

    @Override
    public String toString() {
        return getPfamId();
    }
}

package de.bioforscher.jstructure.model.identifier;

/**
 * The enzyme classification representation.
 * Created by S on 14.07.2017.
 */
public class EnzymeClassificationIdentifier extends AbstractIdentifier {
    private final String ecNumber;

    EnzymeClassificationIdentifier(String ecNumber) {
        this.ecNumber = ecNumber;
    }

    public String getEcNumber() {
        return ecNumber;
    }

    @Override
    public String toString() {
        return getEcNumber();
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;

        EnzymeClassificationIdentifier that = (EnzymeClassificationIdentifier) other;

        return ecNumber.equals(that.ecNumber);
    }

    @Override
    public int hashCode() {
        return ecNumber.hashCode();
    }
}

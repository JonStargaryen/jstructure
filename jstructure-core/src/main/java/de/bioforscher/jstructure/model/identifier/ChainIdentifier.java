package de.bioforscher.jstructure.model.identifier;

/**
 * Represents the name/id of a protein chain. Instances of this class are immutable.
 * Created by bittrich on 4/27/17.
 */
public class ChainIdentifier extends AbstractIdentifier {
    private final ProteinIdentifier proteinIdentifier;
    private final String chainId;
    public static final ChainIdentifier UNKNOWN_CHAIN_IDENTIFIER = IdentifierFactory.createChainIdentifier(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER, "X");

    ChainIdentifier(ProteinIdentifier proteinIdentifier, String chainId) {
        this.proteinIdentifier = proteinIdentifier;
        this.chainId = chainId;
    }

    /**
     * The identifier of the parent of this container.
     * @return e.g. '1brr'
     */
    public ProteinIdentifier getProteinIdentifier() {
        return proteinIdentifier;
    }

    /**
     * The (usually) one-letter-name of this chain.
     * @return e.g. 'A' for the first chain of most proteins
     */
    public String getChainId() {
        return chainId;
    }

    /**
     * The full name of this chain, i.e. pdbId and chainId separated by an underscore.
     * @return e.g. '1brr_A' for the first chain of 1brr
     */
    public String getFullName() {
        return proteinIdentifier.getFullName() + "_" + chainId;
    }

    @Override
    public String toString() {
        return getFullName();
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;

        ChainIdentifier that = (ChainIdentifier) other;

        if (proteinIdentifier != null ? !proteinIdentifier.equals(that.proteinIdentifier) : that.proteinIdentifier != null) return false;
        return chainId != null ? chainId.equals(that.chainId) : that.chainId == null;
    }

    @Override
    public int hashCode() {
        int result = proteinIdentifier != null ? proteinIdentifier.hashCode() : 0;
        result = 31 * result + (chainId != null ? chainId.hashCode() : 0);
        return result;
    }
}

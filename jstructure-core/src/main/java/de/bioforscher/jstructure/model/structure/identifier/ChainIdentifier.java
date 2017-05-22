package de.bioforscher.jstructure.model.structure.identifier;

/**
 * Represents the name/id of a protein chain.
 * Created by bittrich on 4/27/17.
 */
public class ChainIdentifier {
    private ProteinIdentifier pdbId;
    private String chainId;

    public static final ChainIdentifier UNKNOWN_CHAIN_ID = ChainIdentifier.createFromChainId(ProteinIdentifier.UNKNOWN_PROTEIN_ID, "X");

    private ChainIdentifier(ProteinIdentifier pdbId, String chainId) {
        this.pdbId = pdbId;
        this.chainId = chainId;
    }

    /**
     * Creates a new instance of a chainIdentifer.
     * @param pdbId the parent identifier
     * @param chainId this chain's id
     * @return the create instance
     */
    public static ChainIdentifier createFromChainId(ProteinIdentifier pdbId, String chainId) {
        return new ChainIdentifier(pdbId, chainId);
    }

    /**
     * The identifier of the parent of this container.
     * @return e.g. '1brr'
     */
    public ProteinIdentifier getPdbId() {
        return pdbId;
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
        return pdbId.getFullName() + "_" + chainId;
    }

    @Override
    public String toString() {
        return getFullName();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ChainIdentifier that = (ChainIdentifier) o;

        if (pdbId != null ? !pdbId.equals(that.pdbId) : that.pdbId != null) return false;
        return chainId != null ? chainId.equals(that.chainId) : that.chainId == null;
    }

    @Override
    public int hashCode() {
        int result = pdbId != null ? pdbId.hashCode() : 0;
        result = 31 * result + (chainId != null ? chainId.hashCode() : 0);
        return result;
    }
}

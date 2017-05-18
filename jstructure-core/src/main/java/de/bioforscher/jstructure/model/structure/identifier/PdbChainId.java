package de.bioforscher.jstructure.model.structure.identifier;

/**
 * Represents the name/id of a protein chain.
 * Created by bittrich on 4/27/17.
 */
public class PdbChainId {
    //TODO pattern for chainIds?
    private PdbId pdbId;
    private String chainId;

    public static final PdbChainId UNKNOWN_CHAIN_ID = PdbChainId.createFromChainId(PdbId.UNKNOWN_PROTEIN_ID, "X");

    private PdbChainId(PdbId pdbId, String chainId) {
        this.pdbId = pdbId;
        this.chainId = chainId;
    }

    /**
     *
     * @param pdbId
     * @param chainId
     * @return
     */
    public static PdbChainId createFromChainId(PdbId pdbId, String chainId) {
        return new PdbChainId(pdbId, chainId);
    }

    /**
     *
     * @return
     */
    public PdbId getPdbId() {
        return pdbId;
    }

    /**
     *
     * @return
     */
    public String getChainId() {
        return chainId;
    }

    /**
     *
     * @return
     */
    public String getFullName() {
        return pdbId.getFullName() + "_" + chainId;
    }

    @Override
    public String toString() {
        return getFullName();
    }
}

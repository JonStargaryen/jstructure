package de.bioforscher.jstructure.feature.sse;

/**
 * Container that represents a beta Bridge between two residues. It contains the
 * two partner indices and the type of the bridge. For consistency, partner1 is
 * always the small index.
 *
 * @author Aleix Lafita
 *
 */
public class BetaBridge {
    private BridgeType type;
    private int partner1;
    private int partner2;

    public BetaBridge(int i, int j, BridgeType t) {
        this.partner1 = Math.min(i, j);
        this.partner2 = Math.max(i, j);
        this.type = t;
    }

    public BridgeType getType() {
        return this.type;
    }

    public int getPartner1() {
        return this.partner1;
    }

    public int getPartner2() {
        return this.partner2;
    }
}
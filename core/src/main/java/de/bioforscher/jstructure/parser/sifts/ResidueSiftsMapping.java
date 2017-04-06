package de.bioforscher.jstructure.parser.sifts;

/**
 * Residue-specfic SIFTS-mappings.
 * Created by bittrich on 4/6/17.
 */
public class ResidueSiftsMapping {
    private String uniProtId;
    private int uniProtResidueNumber;

    ResidueSiftsMapping(String uniProtId, String uniProtResidueNumber) {
        this.uniProtId = uniProtId;
        this.uniProtResidueNumber = Integer.valueOf(uniProtResidueNumber);
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public int getUniProtResidueNumber() {
        return uniProtResidueNumber;
    }
}

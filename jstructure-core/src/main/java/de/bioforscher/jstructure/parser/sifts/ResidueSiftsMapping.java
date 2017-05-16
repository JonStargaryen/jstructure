package de.bioforscher.jstructure.parser.sifts;

/**
 * Residue-specfic SIFTS-mappings.
 * Created by bittrich on 4/6/17.
 */
public class ResidueSiftsMapping {
    private String uniProtId;
    private int uniProtResidueNumber;

    public static final ResidueSiftsMapping MISSING_VALUE = new ResidueSiftsMapping("?", Integer.MIN_VALUE);

    ResidueSiftsMapping(String uniProtId, String uniProtResidueNumber) {
        this(uniProtId, Integer.valueOf(uniProtResidueNumber));
    }

    ResidueSiftsMapping(String uniProtId, int uniProtResidueNumber) {
        this.uniProtId = uniProtId;
        this.uniProtResidueNumber = uniProtResidueNumber;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public int getUniProtResidueNumber() {
        return uniProtResidueNumber;
    }
}

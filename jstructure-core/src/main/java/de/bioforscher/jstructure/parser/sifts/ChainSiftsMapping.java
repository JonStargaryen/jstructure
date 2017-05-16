package de.bioforscher.jstructure.parser.sifts;

/**
 * Chain-specific SIFTS-mappings.
 * Created by bittrich on 4/6/17.
 */
public class ChainSiftsMapping {
    private String uniProtId, ecNumber, pfam;

    public ChainSiftsMapping(String uniProtId, String ecNumber, String pfam) {
        this.uniProtId = uniProtId;
        this.ecNumber = ecNumber;
        this.pfam = pfam;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public String getEcNumber() {
        return ecNumber;
    }

    public String getPfam() {
        return pfam;
    }
}

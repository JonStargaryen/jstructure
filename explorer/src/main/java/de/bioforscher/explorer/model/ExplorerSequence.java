package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.sifts.ChainSiftsMapping;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;

/**
 * An aligned sequence.
 * Created by bittrich on 4/6/17.
 */
public class ExplorerSequence {
    private String id, title, ec, pfam, uniprot, sequence;

    public ExplorerSequence() {
    }

    public ExplorerSequence(Chain chain) {
        this.id = chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId();
        this.title = chain.getParentProtein().getTitle();
        ChainSiftsMapping chainSiftsMapping = chain.getFeature(ChainSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING);
        this.ec = chainSiftsMapping.getEcNumber();
        this.pfam = chainSiftsMapping.getPfam();
        this.uniprot = chainSiftsMapping.getUniProtId();
        this.sequence = chain.getAminoAcidSequence();
    }

    public String getId() {
        return id;
    }

    public String getTitle() {
        return title;
    }

    public String getEc() {
        return ec;
    }

    public String getPfam() {
        return pfam;
    }

    public String getUniprot() {
        return uniprot;
    }

    public String getSequence() {
        return sequence;
    }
}

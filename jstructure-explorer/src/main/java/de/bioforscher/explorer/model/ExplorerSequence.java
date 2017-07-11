package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.mapping.ChainMapping;
import de.bioforscher.jstructure.model.structure.Chain;

/**
 * An aligned sequence.
 * Created by bittrich on 4/6/17.
 */
@Deprecated
public class ExplorerSequence {
    private String id, title, ec, pfam, uniprot, sequence;

    public ExplorerSequence() {
    }

    public ExplorerSequence(Chain chain) {
        this.id = chain.getChainIdentifier().getFullName();
        this.title = chain.getParentProtein().getTitle();
        this.sequence = chain.getAminoAcidSequence();
        try {
            ChainMapping chainSiftsMapping = chain.getFeatureContainer().getFeature(ChainMapping.class);
            this.ec = chainSiftsMapping.getEcNumber();
            this.pfam = chainSiftsMapping.getPfamId();
            this.uniprot = chainSiftsMapping.getUniProtId();
        } catch (NullPointerException e) {

        }
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

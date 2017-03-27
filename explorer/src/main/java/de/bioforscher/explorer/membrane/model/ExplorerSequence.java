package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;

/**
 * An aligned sequence.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerSequence {
    private String id, title, ec, pfam, uniprot;

    public ExplorerSequence() {
    }

    public ExplorerSequence(Chain chain) {
        this.id = ExplorerModelFactory.getGlobalId(chain);
        this.title = chain.getParentProtein().getTitle();
        this.ec = chain.getFeature(String.class, SiftsParser.EC_NUMBER);
        this.pfam = chain.getFeature(String.class, SiftsParser.PFAM_ID);
        this.uniprot = chain.getFeature(String.class, SiftsParser.UNIPROT_ID);
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
}

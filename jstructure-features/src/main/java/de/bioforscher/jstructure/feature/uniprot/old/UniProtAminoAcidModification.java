package de.bioforscher.jstructure.feature.uniprot.old;

import org.jsoup.nodes.Element;

/**
 * Annotation of modified amino acids.
 * Created by bittrich on 4/18/17.
 */
@Deprecated
public class UniProtAminoAcidModification {
    private int position;
    private String description;
    private Type type;

    public enum Type {
        GLYCOSYLATION_SITE,
        LIPID_BINDING,
        MODIFIED_RESIDUE
    }

    public UniProtAminoAcidModification() {
    }

    UniProtAminoAcidModification(Element describingElement) {
        this.position = Integer.valueOf(describingElement.getElementsByTag("position").first().attr("position"));
        this.description = describingElement.attr("description");
        this.type = mapType(describingElement.attr("type"));
    }

    private Type mapType(String type) {
        if(type.equals("glycosylation site")) {
            return Type.GLYCOSYLATION_SITE;
        }
        if(type.equals("lipid moiety-binding region")) {
            return Type.LIPID_BINDING;
        }
        if(type.equals("modified residue")) {
            return Type.MODIFIED_RESIDUE;
        }

        throw new IllegalArgumentException("unexpected modification type: " + type);
    }

    public int getPosition() {
        return position;
    }

    public String getDescription() {
        return description;
    }

    public Type getType() {
        return type;
    }
}

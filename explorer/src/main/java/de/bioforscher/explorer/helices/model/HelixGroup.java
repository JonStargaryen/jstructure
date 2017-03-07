package de.bioforscher.explorer.helices.model;

import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 3/7/17.
 */
public class HelixGroup {
    private int resn;
    private List<HelixAtom> atoms;
    private String olc, tlc, type;
    private double rasa, sse;

    public HelixGroup() {

    }

    HelixGroup(Group group) {
        this.resn = group.getResidueNumber();
        this.atoms = group.atoms()
                .map(HelixAtom::new)
                .collect(Collectors.toList());
        this.olc = group.getGroupInformation().getOneLetterCode();
        this.tlc = group.getGroupInformation().getThreeLetterCode();

        if(group.isAminoAcid()) {
            this.type = "aa";
        } else if(group.isNucleotide()) {
            this.type = "nucl";
        } else {
            this.type = "het";
        }
    }

    public int getResn() {
        return resn;
    }

    public List<HelixAtom> getAtoms() {
        return atoms;
    }

    public String getOlc() {
        return olc;
    }

    public String getTlc() {
        return tlc;
    }

    public String getType() {
        return type;
    }

    public double getSse() {
        return sse;
    }

    public double getRasa() {
        return rasa;
    }
}

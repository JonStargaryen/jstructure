package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * An aligned sequence.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerSequence {
    private String id;
    private List<ExplorerAminoAcid> sequence;

    public ExplorerSequence() {
    }

    public ExplorerSequence(Chain chain, Map<String, String> alignedSequences) {
        this.id = ExplorerModelFactory.getGlobalId(chain);

        this.sequence = new ArrayList<>();
        String alignedSequence = alignedSequences.get(id);
        Iterator<Group> aminoAcids = chain.aminoAcids().iterator();
        for(int pos = 0; pos < alignedSequence.length(); pos++) {
            char olc = alignedSequence.charAt(pos);
            String olcString = String.valueOf(olc);
            ExplorerAminoAcid aminoAcid;
            if(olc != '-') {
                aminoAcid = new ExplorerAminoAcid(pos, olcString, aminoAcids.next());
            } else {
                aminoAcid = new ExplorerAminoAcid(pos, olcString);
            }
            this.sequence.add(aminoAcid);
        }
    }

    public String getId() {
        return id;
    }

    public List<ExplorerAminoAcid> getSequence() {
        return sequence;
    }
}

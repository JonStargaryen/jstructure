package de.bioforscher.explorer.membrane.model;

import java.util.ArrayList;
import java.util.List;

/**
 * An aligned sequence.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerSequence {
    private String id;
    private List<SequencePosition> sequence;

    public ExplorerSequence() {
    }

    public ExplorerSequence(ExplorerChain chain, String id, String sequence) {
        this.id = id;

        this.sequence = new ArrayList<>();
        int groupIndex = 0;
        for(int pos = 0; pos < sequence.length(); pos++) {
            String olc = sequence.substring(pos, pos + 1);
            this.sequence.add(new SequencePosition(chain.getAminoAcids().get(groupIndex).getResn(), olc));
            // increase 'pointer' when non-gap position was observed
            if(!olc.equals("-")) {
                groupIndex++;
            }
        }
    }

    public String getId() {
        return id;
    }

    public List<SequencePosition> getSequence() {
        return sequence;
    }
}

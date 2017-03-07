package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

/**
 * The representation of helices and kinks therein.
 * Created by bittrich on 2/28/17.
 */
public class ExplorerHelix {
    private static final int GRACE_INTERVAL = 1;
    private String chain;
    private int start, end;

    public ExplorerHelix() {

    }

    public ExplorerHelix(String chain, int start, int end) {
        this.chain = chain;
        this.start = start;
        this.end = end;
    }


    public String getChain() {
        return chain;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    static Stream<ExplorerHelix> extract(Chain chain) {
//        List<Integer> helixPositions = chain.groups()
//                .filter(group -> group.getFeature(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure().isHelixType())
//                .map(Group::getResidueNumber)
//                .collect(Collectors.toList());

        List<ExplorerHelix> helices = new ArrayList<>();

        int helixStart = -1;
        int grace = 0;
        for(Group group : chain.getGroups()) {
            if(!group.isAminoAcid()) {
                continue;
            }
            boolean isHelix = isHelix(group);

            if(helixStart == -1 && isHelix) {
                // helix start
                helixStart = group.getResidueNumber();
                continue;
            }

            if(isHelix) {
                // helix is extended
                grace = 0;
                continue;
            }

            if(helixStart != -1 && !isHelix) {
                // helix end?

                //TODO detect actually turning helices and separate them
                if(grace < GRACE_INTERVAL) {
                    // 1 non-helix residue grace to connect helices
                    grace++;
                    continue;
                }

                helices.add(new ExplorerHelix(chain.getChainId(), helixStart, group.getResidueNumber() - (grace + 1)));
                helixStart = -1;
                grace = 0;
            }
        }

        return helices.stream();
    }

    private static boolean isHelix(Group group) {
        return group.getFeature(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure().isHelixType();
    }
}

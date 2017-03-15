package de.bioforscher.explorer.membrane.model.homologous;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;

import java.util.stream.Collectors;

/**
 * Neutered chain.
 * Created by bittrich on 3/15/17.
 */
public class HomologousChain {
    private String id;
    private String sequence;

    public HomologousChain() {
    }

    public HomologousChain(Chain chain) {
        this.id = chain.getChainId();
        this.sequence = chain.aminoAcids()
                .map(Group::getGroupInformation)
                .map(GroupInformation::getOneLetterCode)
                .collect(Collectors.joining());
    }

    public String getId() {
        return id;
    }

    public String getSequence() {
        return sequence;
    }
}

package de.bioforscher.explorer.helices.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 3/7/17.
 */
public class HelixChain {
    private String id;
    private List<HelixGroup> groups;

    public HelixChain() {

    }

    HelixChain(Chain chain) {
        this.id = chain.getChainId();
        this.groups = chain.groups()
                .map(HelixGroup::new)
                .collect(Collectors.toList());
    }

    public String getId() {
        return id;
    }

    public List<HelixGroup> getGroups() {
        return groups;
    }
}

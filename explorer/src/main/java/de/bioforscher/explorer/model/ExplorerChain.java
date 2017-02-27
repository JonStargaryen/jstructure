package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The reduced representation of {@link Chain} objects.
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerChain {
    private String id;
    private List<ExplorerGroup> groups;

    public ExplorerChain() {

    }

    ExplorerChain(Chain chain) {
        this.id = chain.getChainId();
        this.groups = chain.groups()
                .map(ExplorerGroup::new)
                .collect(Collectors.toList());
    }

    public String getId() {
        return id;
    }

    public List<ExplorerGroup> getGroups() {
        return groups;
    }
}

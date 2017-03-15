package de.bioforscher.explorer.membrane.model.homologous;

import de.bioforscher.jstructure.model.structure.Protein;

import java.util.List;
import java.util.stream.Collectors;

/**
 * One concrete observation of an homologous protein.
 * Created by bittrich on 3/15/17.
 */
public class HomologousProtein {
    private String name, title;
    private List<HomologousChain> chains;

    public HomologousProtein() {
    }

    public HomologousProtein(Protein protein) {
        this.name = protein.getName();
        this.title = protein.getTitle();
        this.chains = protein.chains()
                .filter(chain -> chain.aminoAcids().count() > 0)
                .map(HomologousChain::new)
                .collect(Collectors.toList());
    }

    public String getName() {
        return name;
    }

    public String getTitle() {
        return title;
    }

    public List<HomologousChain> getChains() {
        return chains;
    }
}

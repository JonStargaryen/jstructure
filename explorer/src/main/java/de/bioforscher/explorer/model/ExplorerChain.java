package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.uniprot.*;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The reduced representation of {@link Chain} objects.
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerChain {
    private List<UniProtReference> references;
    private List<UniProtNaturalVariant> variants;
    private List<UniProtMutagenesisSite> mutations;
    private String id;
    private List<ExplorerGroup> groups;

    public ExplorerChain() {

    }

    ExplorerChain(Chain chain) {
        this.id = chain.getChainId();
        this.groups = chain.groups()
                .map(ExplorerGroup::new)
                .collect(Collectors.toList());

        UniProtAnnotationContainer uniProtAnnotationContainer = chain.getFeature(UniProtAnnotationContainer.class, UniProtAnnotator.UNIPROT_ANNOTATION);
        this.mutations = uniProtAnnotationContainer.getMutagenesisSites();
        this.variants = uniProtAnnotationContainer.getNaturalVariants();
        this.references = uniProtAnnotationContainer.getReferences();
    }

    public String getId() {
        return id;
    }

    public List<ExplorerGroup> getGroups() {
        return groups;
    }

    public List<UniProtReference> getReferences() {
        return references;
    }

    public List<UniProtNaturalVariant> getVariants() {
        return variants;
    }

    public List<UniProtMutagenesisSite> getMutations() {
        return mutations;
    }
}

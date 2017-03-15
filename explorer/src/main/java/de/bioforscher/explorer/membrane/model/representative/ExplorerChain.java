package de.bioforscher.explorer.membrane.model.representative;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
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
    private String id, uniprot, pfam, ec;
    private List<ExplorerGroup> groups;
    private boolean hasAminoAcids;

    public ExplorerChain() {

    }

    ExplorerChain(Chain chain) {
        this.id = chain.getChainId();
        this.groups = chain.groups()
                .map(ExplorerGroup::new)
                .collect(Collectors.toList());

        this.uniprot = chain.getFeature(String.class, SiftsParser.UNIPROT_ID);
        this.pfam = chain.getFeature(String.class, SiftsParser.PFAM_ID);
        this.ec = chain.getFeature(String.class, SiftsParser.EC_NUMBER);

        try {
            UniProtAnnotationContainer uniProtAnnotationContainer = chain.getFeature(UniProtAnnotationContainer.class, UniProtAnnotator.UNIPROT_ANNOTATION);
            this.mutations = uniProtAnnotationContainer.getMutagenesisSites();
            this.variants = uniProtAnnotationContainer.getNaturalVariants();
            this.references = uniProtAnnotationContainer.getReferences();
        } catch (NullPointerException e) {
            //TODO unify this pattern - featureMap returning optionals?
        }

        this.hasAminoAcids = Selection.on(chain)
                .aminoAcids()
                .asFilteredGroups()
                .count() > 0;
    }

    public String getId() {
        return id;
    }

    public String getUniprot() {
        return uniprot;
    }

    public String getPfam() {
        return pfam;
    }

    public String getEc() {
        return ec;
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

    public boolean isHasAminoAcids() {
        return hasAminoAcids;
    }
}

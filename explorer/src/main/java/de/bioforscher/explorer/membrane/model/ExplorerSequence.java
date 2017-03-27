package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import de.bioforscher.jstructure.parser.uniprot.*;

import java.util.*;

/**
 * An aligned sequence.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerSequence {
    private String id, title, ec, pfam, uniprot;
    private List<ExplorerAminoAcid> sequence;
    private List<UniProtMutagenesisSite> mutations;
    private List<UniProtNaturalVariant> variants;
    private List<UniProtReference> references;

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

        this.title = chain.getParentProtein().getTitle();
        this.ec = chain.getFeature(String.class, SiftsParser.EC_NUMBER);
        this.pfam = chain.getFeature(String.class, SiftsParser.PFAM_ID);
        this.uniprot = chain.getFeature(String.class, SiftsParser.UNIPROT_ID);

        this.references = new ArrayList<>();
        try {
            UniProtAnnotationContainer uniProtAnnotationContainer = chain.getFeature(UniProtAnnotationContainer.class, UniProtAnnotator.UNIPROT_ANNOTATION);
            this.mutations = uniProtAnnotationContainer.getMutagenesisSites();
            this.variants = uniProtAnnotationContainer.getNaturalVariants();
            mutations.stream()
                    .map(UniProtMutagenesisSite::getEvidence)
                    .flatMap(Collection::stream)
                    .map(evidenceNumber -> uniProtAnnotationContainer.getReferences().get(Integer.valueOf(evidenceNumber) - 1))
                    .forEach(references::add);
            variants.stream()
                    .map(UniProtNaturalVariant::getEvidence)
                    .flatMap(Collection::stream)
                    .map(evidenceNumber -> uniProtAnnotationContainer.getReferences().get(Integer.valueOf(evidenceNumber) - 1))
                    .forEach(references::add);
        } catch (NullPointerException e) {
            this.mutations = new ArrayList<>();
            this.variants = new ArrayList<>();
        }
    }

    public String getId() {
        return id;
    }

    public List<ExplorerAminoAcid> getSequence() {
        return sequence;
    }

    public String getTitle() {
        return title;
    }

    public String getEc() {
        return ec;
    }

    public String getPfam() {
        return pfam;
    }

    public String getUniprot() {
        return uniprot;
    }

    public List<UniProtMutagenesisSite> getMutations() {
        return mutations;
    }

    public List<UniProtNaturalVariant> getVariants() {
        return variants;
    }

    public List<UniProtReference> getReferences() {
        return references;
    }
}

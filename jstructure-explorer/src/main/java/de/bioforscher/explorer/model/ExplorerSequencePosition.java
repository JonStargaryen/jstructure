package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.uniprot.UniProtActiveSite;
import de.bioforscher.jstructure.feature.uniprot.UniProtAnnotationContainer;
import de.bioforscher.jstructure.feature.uniprot.UniProtMutagenesisSite;
import de.bioforscher.jstructure.feature.uniprot.UniProtNaturalVariant;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * One position of the MSA-view.
 * Created by bittrich on 4/6/17.
 */
@Deprecated
public class ExplorerSequencePosition {
    private String olcs;
    private boolean uniform, mutant, variant, activeSite;
    private List<UniProtMutagenesisSite> mutations;
    private List<UniProtNaturalVariant> variants;
    private List<UniProtActiveSite> activeSites;

    public ExplorerSequencePosition() {
    }

    public ExplorerSequencePosition(Set<UniProtAnnotationContainer> containers, List<Chain> chains, int position) {
        this.olcs = chains.stream()
                .map(chain -> chain.select().aminoAcids().residueNumber(position).asOptionalGroup())
                .map(optionalGroup -> optionalGroup.map(group -> ((AminoAcid) group).getOneLetterCode()).orElse("-"))
                .collect(Collectors.joining(System.lineSeparator()));

        // determine if variant positions
        this.uniform = Pattern.compile(System.lineSeparator()).splitAsStream(olcs)
                .distinct()
                .count() == 1;

        String positionString = String.valueOf(position);
        this.mutations = containers.stream()
                .map(UniProtAnnotationContainer::getMutagenesisSites)
                .flatMap(Collection::stream)
                .distinct()
                .filter(mutagenesisSite -> mutagenesisSite.getPosition().contains(positionString))
                .collect(Collectors.toList());
        this.variants = containers.stream()
                .map(UniProtAnnotationContainer::getNaturalVariants)
                .flatMap(Collection::stream)
                .distinct()
                .filter(naturalVariant -> naturalVariant.getPosition().contains(positionString))
                .collect(Collectors.toList());
        this.activeSites = containers.stream()
                .map(UniProtAnnotationContainer::getActiveSites)
                .flatMap(Collection::stream)
                .distinct()
                .filter(activeSite -> activeSite.getPosition().equals(positionString))
                .collect(Collectors.toList());
        this.mutant = !mutations.isEmpty();
        this.variant = !variants.isEmpty();
        this.activeSite = !activeSites.isEmpty();
    }

    public String getOlcs() {
        return olcs;
    }

    public boolean isUniform() {
        return uniform;
    }

    public boolean isMutant() {
        return mutant;
    }

    public boolean isVariant() {
        return variant;
    }

    public boolean isActiveSite() {
        return activeSite;
    }

    public List<UniProtMutagenesisSite> getMutations() {
        return mutations;
    }

    public List<UniProtNaturalVariant> getVariants() {
        return variants;
    }

    public List<UniProtActiveSite> getActiveSites() {
        return activeSites;
    }
}

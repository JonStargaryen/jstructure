package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotationContainer;
import de.bioforscher.jstructure.parser.uniprot.UniProtMutagenesisSite;
import de.bioforscher.jstructure.parser.uniprot.UniProtNaturalVariant;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * One position of the MSA.
 * Created by bittrich on 3/27/17.
 */
public class ExplorerSequencePosition {
    private String olcs;
    private boolean uniform, mutant, variant;
    private List<UniProtMutagenesisSite> mutations;
    private List<UniProtNaturalVariant> variants;

    public ExplorerSequencePosition() {
    }

    public ExplorerSequencePosition(Set<UniProtAnnotationContainer> containers, Map<String, String> alignedSequences, int position) {
        this.olcs = alignedSequences.values().stream()
                .map(sequence -> sequence.charAt(position))
                .map(String::valueOf)
                .collect(Collectors.joining(System.lineSeparator()));

        // determine if variant positions
        this.uniform = Pattern.compile(System.lineSeparator()).splitAsStream(olcs)
                .distinct()
                .count() == 1;

        String positionString = String.valueOf(position);
        //TODO there is a offset between this position (renumbered) and the uniProt data
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
        this.mutant = !mutations.isEmpty();
        this.variant = !variants.isEmpty();
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

    public List<UniProtMutagenesisSite> getMutations() {
        return mutations;
    }

    public List<UniProtNaturalVariant> getVariants() {
        return variants;
    }
}

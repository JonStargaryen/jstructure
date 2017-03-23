package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represent the multi-sequence alignment of 1 cluster.
 * Created by bittrich on 3/17/17.
 */
public class ExplorerAlignment {
    private String representativeChainId;
    private List<String> homologous;
    private List<ExplorerSequence> sequences;
    private List<Integer> variants;

    public ExplorerAlignment() {
    }

    public ExplorerAlignment(String representativeChainId, List<Chain> chains, Map<String, String> alignedSequences) {
        this.representativeChainId = representativeChainId;
        this.homologous = alignedSequences.keySet().stream().collect(Collectors.toList());
        this.sequences = chains.stream()
                .map(chain -> new ExplorerSequence(chain, alignedSequences))
                .collect(Collectors.toList());
        // determine variant positions
        this.variants = IntStream.range(0, sequences.get(0).getSequence().size())
                .filter(pos -> alignedSequences
                        .values()
                        .stream()
                        .map(sequence -> sequence.charAt(pos))
                        .distinct()
                        .count() > 1)
                .boxed()
                .collect(Collectors.toList());
    }

    public String getRepresentativeChainId() {
        return representativeChainId;
    }

    public List<ExplorerSequence> getSequences() {
        return sequences;
    }

    public List<String> getHomologous() {
        return homologous;
    }

    public List<Integer> getVariants() {
        return variants;
    }
}
